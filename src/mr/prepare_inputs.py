"""Module 1 data step: derive analysis-ready, cis-region subsets of the raw
downloads (src/mr/download_inputs.sh) for the exposure and the positive-control
traits. Nothing here touches an Alzheimer's outcome file.

    uv run --frozen python -m src.mr.prepare_inputs

Why a derived file at all:
- eQTLGen distributes Z-scores and the assessed-allele frequency in two
  separate files; src/mr/sumstats.R needs both in one table. The effect-allele
  frequency (eaf) is joined by rsID and allele check here.
- data.table::fread needs the R.utils package to read .gz, which is not in the
  renv library; the subsets are written uncompressed.
- Some control GWAS have no usable rsID column (Evangelou 2018 BP; MAGIC 2021
  fasting glucose; the resting-heart-rate file mixes rsIDs and chr:pos ids).
  rsIDs are attached by matching chromosome, GRCh37 position AND the unordered
  allele pair to the 1000 Genomes EUR LD panel's .bim (the ID space clumping
  uses). A position with no matching panel variant, or with several, is
  dropped and counted.

Rows are restricted to each target gene's Ensembl GRCh37 body +/- the
pre-registered sensitivity window (window_kb_sensitivity, 1000 kb), plus the
range of that gene's eQTLGen cis rows. That contains every possible
instrument, so the subset is a speed measure, not an analysis filter: the
pipeline applies its own +/-100 kb (primary) or +/-1000 kb (sensitivity) window
afterwards. The original column names are kept (new ones are appended), so
config/mr_targets.yaml maps to the published headers. Every row count (read,
in region, rsID mapped / unmapped / ambiguous) is logged.

Needs: data/interim/gene_regions_GRCh37.csv (src/mr/regions.R) and PLINK 1.9
on PATH (allele frequencies of the 1000G EUR panel).
"""

from __future__ import annotations

import csv
import io
import subprocess
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import yaml

from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.common.runlog import start_run

RAW = REPO_ROOT / "data" / "raw" / "mr"
OUT = REPO_ROOT / "data" / "interim" / "mr"

EQTL_FILE = RAW / "eqtlgen" / "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
AF_FILE = RAW / "eqtlgen" / "2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
# HMGCR is not a target: it is used only for the informal LDL -> CAD pipeline
# check (the formal validation target needs researcher-transcribed inputs).
EXTRA_GENES = ["HMGCR"]

# trait -> source file, kind of reader, and the columns that hold chrom, pos,
# allele 1 and allele 2 (for region filter and rsID mapping) and rsID (None
# = map from the 1000G EUR panel).
TRAITS = {
    # Header: SNP CHR BP A1 A0 MAF HWEP INFO BETA SE P (readme: A1 = effect allele).
    # Only the MINOR-allele frequency is given, so the A1 frequency is oriented
    # with the 1000G EUR frequency of that allele (non-palindromic variants only:
    # for A/T and C/G variants MAF cannot say which allele is minor without
    # assuming the strand, so their eaf is left NA and the reader drops them).
    # SNP ids are a mix of rsIDs and chr:pos_A_B, so rsIDs come from the panel.
    "resting_heart_rate": {
        "path": RAW / "gwascatalog" / "GCST007609" / "ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz",
        "reader": "gz", "chrom": "CHR", "pos": "BP", "a1": "A1", "a2": "A0", "rsid": None, "maf": "MAF"},
    "diastolic_blood_pressure": {
        "path": RAW / "gwascatalog" / "GCST006630" / "Evangelou_30224653_DBP.txt.gz",
        "reader": "gz", "chrom": None, "pos": None, "markername": "MarkerName", "a1": "Allele1", "a2": "Allele2", "rsid": None},
    "systolic_blood_pressure": {
        "path": RAW / "gwascatalog" / "GCST006624" / "Evangelou_30224653_SBP.txt.gz",
        "reader": "gz", "chrom": None, "pos": None, "markername": "MarkerName", "a1": "Allele1", "a2": "Allele2", "rsid": None},
    # effect_allele_frequency (1000G phase 1 v3, per the README) is NA for some
    # rows; eaf_filled keeps it where present and otherwise takes the 1000G
    # Phase 3 EUR frequency of the effect allele - non-palindromic variants
    # only (palindromic ones stay NA and are dropped by the reader).
    "fasting_glucose": {
        "path": RAW / "gwascatalog" / "GCST90002232" / "GCST90002232_buildGRCh37.tsv.gz",
        "reader": "gz", "chrom": "chromosome", "pos": "base_pair_location", "a1": "effect_allele", "a2": "other_allele",
        "rsid": None, "eaf_fill": "effect_allele_frequency"},
    "sleep_duration": {
        "path": RAW / "gwascatalog" / "GCST007561" / "Dashti_30846698_sleepdurationsumstats.txt.zip",
        "reader": "zip", "chrom": "CHR", "pos": "BP", "a1": "ALLELE1", "a2": "ALLELE0", "rsid": "SNP"},
    "ldl_cholesterol": {
        "path": RAW / "gwascatalog" / "GCST002222" / "jointGwasMc_LDL.txt.gz",
        "reader": "gz", "chrom": None, "pos": None, "chrpos": "SNP_hg19", "a1": "A1", "a2": "A2", "rsid": "rsid"},
    "coronary_artery_disease": {
        "path": RAW / "gwascatalog" / "GCST003116" / "cad.add.160614.website.txt",
        "reader": "plain", "chrom": "chr", "pos": "bp_hg19", "a1": "effect_allele", "a2": "noneffect_allele", "rsid": "markername"},
}


def close_stream(proc) -> None:
    if proc is not None and proc.wait() != 0:
        raise RuntimeError(f"decompression failed (exit {proc.returncode}): file truncated or corrupt?")


def open_text(path: Path, reader: str):
    if reader == "gz":
        p = subprocess.Popen(["zcat", str(path)], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        return p.stdout, p
    if reader == "zip":
        z = zipfile.ZipFile(path)
        members = [m for m in z.namelist() if not m.endswith("/") and not m.startswith("__MACOSX")]
        if len(members) != 1:
            raise RuntimeError(f"{path}: expected one member in the zip, found {members}")
        return io.TextIOWrapper(z.open(members[0]), encoding="utf-8"), None
    return open(path, encoding="utf-8"), None


def sniff_split(header_line: str) -> str:
    return "\t" if "\t" in header_line else None  # None = any whitespace


def target_genes() -> list[str]:
    cfg = yaml.safe_load((CONFIG_DIR / "mr_targets.yaml").read_text(encoding="utf-8"))
    return list(cfg["targets"]) + EXTRA_GENES


def extract_eqtlgen(genes: list[str], run) -> tuple[Path, dict]:
    """eQTLGen rows for the genes, with the assessed-allele frequency joined."""
    out = OUT / "eqtlgen_cis_targets_with_eaf.tsv"
    gene_set = set(genes)
    rows = []
    stream, proc = open_text(EQTL_FILE, "gz")
    header = stream.readline().rstrip("\n").split("\t")
    need = ["Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele", "Zscore", "Gene",
            "GeneSymbol", "NrSamples"]
    miss = [c for c in need if c not in header]
    if miss:
        raise RuntimeError(f"eQTLGen header lacks {miss}: {header}")
    ig = header.index("GeneSymbol")
    n_total = 0
    for line in stream:
        n_total += 1
        f = line.rstrip("\n").split("\t")
        if f[ig] in gene_set:
            rows.append(f)
    close_stream(proc)
    run.info(f"eQTLGen: {n_total} rows read, {len(rows)} for {len(gene_set)} genes")
    snps = {r[header.index("SNP")] for r in rows}

    af = {}
    stream, proc = open_text(AF_FILE, "gz")
    ah = stream.readline().rstrip("\n").split("\t")
    for c in ["SNP", "hg19_chr", "hg19_pos", "AlleleA", "AlleleB", "AlleleB_all"]:
        if c not in ah:
            raise RuntimeError(f"eQTLGen AF header lacks {c}: {ah}")
    i_snp, i_a, i_b, i_f = ah.index("SNP"), ah.index("AlleleA"), ah.index("AlleleB"), ah.index("AlleleB_all")
    for line in stream:
        f = line.rstrip("\n").split("\t")
        if f[i_snp] in snps:
            af[f[i_snp]] = (f[i_a].upper(), f[i_b].upper(), f[i_f], f[ah.index("hg19_chr")], f[ah.index("hg19_pos")])
    close_stream(proc)
    # Build check: eQTLGen SNPPos must equal the AF file's hg19 position.
    i_s, i_c, i_p = header.index("SNP"), header.index("SNPChr"), header.index("SNPPos")
    same = sum(1 for r in rows if r[i_s] in af and (af[r[i_s]][3], af[r[i_s]][4]) == (r[i_c], r[i_p]))
    have = sum(1 for r in rows if r[i_s] in af)
    run.info(f"eQTLGen build check: {same} of {have} rows have SNPChr/SNPPos equal to the AF file's hg19_chr/hg19_pos")

    ix = {c: header.index(c) for c in need}
    counts = defaultdict(int)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header + ["AssessedAllele_freq", "eaf_source"])
        for f in rows:
            snp = f[ix["SNP"]]
            ea, oa = f[ix["AssessedAllele"]].upper(), f[ix["OtherAllele"]].upper()
            a = af.get(snp)
            if a is None:
                eaf, src = "NA", "no AF record"
            elif a[1] == ea and a[0] == oa:
                eaf, src = a[2], "AlleleB = assessed"
            elif a[0] == ea and a[1] == oa:
                eaf, src = (f"{1 - float(a[2]):.6g}" if a[2] not in ("", "NA") else "NA"), "AlleleA = assessed (1 - AF)"
            else:
                eaf, src = "NA", "AF alleles do not match"
            if eaf in ("", "NA"):
                eaf = "NA"
            counts[src] += 1
            f.extend([eaf, src])  # kept on the row: regions_and_index() reads the eaf back
            w.writerow(f)
    for k, v in sorted(counts.items()):
        run.info(f"eQTLGen eaf join: {k}: {v}")
    per_gene = defaultdict(int)
    for f in rows:
        per_gene[f[ig]] += 1
    for g in genes:
        run.info(f"eQTLGen rows for {g}: {per_gene.get(g, 0)}")
    return out, {"header": header, "rows": rows, "ix": ix, "per_gene": per_gene, "eaf_counts": dict(counts)}


def regions(eq: dict, window_kb: int, run) -> dict:
    """Per gene: list of (chrom, lo, hi): the Ensembl GRCh37 gene body +/- window_kb
    (data/interim/gene_regions_GRCh37.csv, written by src/mr/regions.R) and the
    range of the gene's eQTLGen cis rows."""
    reg = defaultdict(list)
    cache = REPO_ROOT / "data" / "interim" / "gene_regions_GRCh37.csv"
    with cache.open() as fh:
        ens = {r["gene"]: r for r in csv.DictReader(fh)}
    for g in target_genes():
        if g not in ens:
            raise RuntimeError(f"{g} is not in {cache}: run mr_gene_region() for it first")
        r = ens[g]
        reg[g].append((r["chrom"], max(1, int(r["start"]) - window_kb * 1000), int(r["end"]) + window_kb * 1000))
    ix = eq["ix"]
    rng = {}
    for f in eq["rows"]:
        g, c, p = f[ix["GeneSymbol"]], f[ix["SNPChr"]], int(f[ix["SNPPos"]])
        rng[g] = (c, min(rng[g][1], p), max(rng[g][2], p)) if g in rng else (c, p, p)
    for g, r in rng.items():
        reg[g].append(r)
    for g, rs in reg.items():
        run.info(f"region {g}: " + "; ".join(f"chr{c}:{lo}-{hi}" for c, lo, hi in rs))
    return reg


def reference_index(reg: dict, run) -> dict:
    """1000 Genomes EUR (the LD panel) variants in the regions: (chrom, pos) ->
    [(rsid, {allele1, allele2}, {allele: EUR frequency})]. This is the rsID space
    clumping uses, and an exposure-independent source of allele frequencies."""
    bfile = RAW / "ld_1kg_eur" / "g1000_eur"
    by_chrom = defaultdict(list)
    for rs in reg.values():
        for c, lo, hi in rs:
            by_chrom[c].append((lo, hi))
    tmp = OUT / "ref_regions"
    with open(f"{tmp}.ranges", "w") as fh:
        for c, ivs in by_chrom.items():
            for i, (lo, hi) in enumerate(ivs):
                fh.write(f"{c}\t{lo}\t{hi}\tr{c}_{i}\n")
    subprocess.run(["plink", "--bfile", str(bfile), "--extract", "range", f"{tmp}.ranges", "--freq",
                    "--out", str(tmp)], check=True, capture_output=True)
    pos = {}
    with open(bfile.with_suffix(".bim")) as fh:
        for line in fh:
            c, snp, _, p, a1, a2 = line.split()
            if c in by_chrom and in_regions(c, int(p), by_chrom):
                pos[snp] = (c, int(p))
    idx = defaultdict(list)
    with open(f"{tmp}.frq") as fh:
        next(fh)
        for line in fh:
            c, snp, a1, a2, maf, _ = line.split()
            if snp not in pos:
                continue
            f1 = float(maf)
            idx[pos[snp]].append((snp, frozenset({a1.upper(), a2.upper()}), {a1.upper(): f1, a2.upper(): 1 - f1}))
    run.info(f"1000G EUR reference: {sum(len(v) for v in idx.values())} variants in the regions")
    return idx


def in_regions(c: str, p: int, by_chrom: dict) -> bool:
    return any(lo <= p <= hi for lo, hi in by_chrom.get(c, ()))


def extract_trait(name: str, spec: dict, regions: dict, pos_index: dict, run) -> tuple[Path, dict]:
    by_chrom = defaultdict(list)
    for rs in regions.values():
        for c, lo, hi in rs:
            by_chrom[c].append((lo, hi))
    stream, proc = open_text(spec["path"], spec["reader"])
    first = stream.readline().rstrip("\n").rstrip("\r")
    sep = sniff_split(first)
    header = first.split(sep)
    keys = ("chrom", "pos", "a1", "a2", "rsid", "markername", "chrpos", "maf", "eaf_fill")
    for k in keys:
        if spec.get(k) and spec[k] not in header:
            raise RuntimeError(f"{spec['path'].name}: column {spec[k]} not in header {header}")
    gi = {k: header.index(spec[k]) for k in keys if spec.get(k)}
    out = OUT / "positive_controls" / f"{name}_cisregions.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    c = defaultdict(int)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        extra = (["chrom_grch37", "pos_grch37"] + (["rsid_mapped"] if spec["rsid"] is None else [])
                 + (["eaf_oriented"] if spec.get("maf") else [])
                 + (["eaf_filled"] if spec.get("eaf_fill") else []))
        w.writerow(header + extra)
        for line in stream:
            c["read"] += 1
            f = line.rstrip("\n").rstrip("\r").split(sep)
            if len(f) != len(header):
                c["malformed (field count differs from header)"] += 1
                continue
            if "markername" in gi:          # Evangelou: "chr:pos:SNP" / "chr:pos:INDEL"
                parts = f[gi["markername"]].split(":")
                chrom, pos = parts[0], parts[1]
            elif "chrpos" in gi:            # GLGC: "chr1:12345"
                chrom, pos = f[gi["chrpos"]].replace("chr", "").split(":")
            else:
                chrom, pos = f[gi["chrom"]], f[gi["pos"]]
            chrom = chrom.replace("chr", "")
            try:
                pos_i = int(float(pos))
            except ValueError:
                c["position not numeric"] += 1
                continue
            if not in_regions(chrom, pos_i, by_chrom):
                continue
            c["in cis regions"] += 1
            row = f + [chrom, str(pos_i)]
            if spec["rsid"] is None:
                a1, a2 = f[gi["a1"]].upper(), f[gi["a2"]].upper()
                alle = frozenset({a1, a2})
                cand = [(r, fr) for r, al, fr in pos_index.get((chrom, pos_i), []) if al == alle]
                hits = sorted({r for r, _ in cand})
                if len(hits) == 1:
                    c["rsID mapped (chr, pos, alleles match 1000G EUR)"] += 1
                    row.append(hits[0])
                    if spec.get("maf"):
                        fr = cand[0][1]
                        maf = float(f[gi["maf"]])
                        palin = {a1, a2} in ({"A", "T"}, {"C", "G"})
                        if palin:
                            c["eaf NA (palindromic: MAF cannot be oriented without assuming the strand); the reader drops these"] += 1
                            row.append("NA")
                        else:
                            c["eaf oriented from MAF with the 1000G EUR frequency of A1"] += 1
                            row.append(f"{maf:.6g}" if fr[a1] <= 0.5 else f"{1 - maf:.6g}")
                    if spec.get("eaf_fill"):
                        given = f[gi["eaf_fill"]]
                        palin = {a1, a2} in ({"A", "T"}, {"C", "G"})
                        if given not in ("", "NA", "nan"):
                            c["eaf as given in the file"] += 1
                            row.append(given)
                        elif palin:
                            c["eaf NA in file, palindromic: left NA (reader drops)"] += 1
                            row.append("NA")
                        else:
                            c["eaf NA in file, filled with the 1000G EUR frequency of the effect allele"] += 1
                            row.append(f"{cand[0][1][a1]:.6g}")
                elif not hits:
                    c["rsID not mapped (no 1000G EUR variant with same chr, pos, alleles) - row dropped"] += 1
                    continue
                else:
                    c["rsID ambiguous (several 1000G EUR rsIDs) - row dropped"] += 1
                    continue
            w.writerow(row)
            c["written"] += 1
    close_stream(proc)
    for k, v in c.items():
        run.info(f"{name}: {k}: {v}")
    return out, dict(c)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    only = set(sys.argv[1:])
    with start_run(module="mr_prepare_inputs", config_path=CONFIG_DIR / "mr_targets.yaml", seed=None,
                   packages=["pyyaml"]) as run:
        cfg = yaml.safe_load((CONFIG_DIR / "mr_targets.yaml").read_text(encoding="utf-8"))
        genes = target_genes()
        eq_path, eq = extract_eqtlgen(genes, run)
        run.record_output(eq_path)
        for g in genes:
            if eq["per_gene"].get(g, 0) == 0:
                run.warning(f"{g}: no eQTLGen cis rows (gene not among the 19,250 genes eQTLGen tested)")
        # widest window any analysis uses (the pre-registered sensitivity window)
        reg = regions(eq, int(cfg["window_kb_sensitivity"]), run)
        ref = reference_index(reg, run)
        for name, spec in TRAITS.items():
            if only and name not in only:
                continue
            if not spec["path"].exists():
                run.warning(f"{name}: raw file {spec['path']} missing - skipped")
                continue
            path, _ = extract_trait(name, spec, reg, ref, run)
            run.record_output(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
