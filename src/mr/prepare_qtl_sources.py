"""Module 1, added instrument sources (decision 2026-10-21, option b): build the
analysis-ready GRCh37 tables for the plasma pQTL (UKB-PPP) and brain eQTL
(MetaBrain) sources from their OPEN distributions. Nothing here touches an
Alzheimer's outcome file.

    uv run --frozen python -m src.mr.prepare_qtl_sources

Access (checked 2026-10-21; details in config/mr_targets.yaml): the full
UKB-PPP files need a Synapse account and the full MetaBrain files a form, so
neither was downloaded. The open versions used:

* UKB-PPP: Sun et al. 2023 Nature, Supplementary Table 9 (sheet "ST9" of
  41586_2023_6592_MOESM3_ESM.xlsx): every significant primary pQTL of the
  European discovery cohort. Rows whose "Assay Target" is a target gene (or
  HMGCR) are kept, cis and trans (the pipeline's cis window removes trans
  rows). GENPOS (hg38) is lifted to GRCh37 with src/mr/liftover.py; the lifted
  position is checked against the hg37 position inside the variant id, and the
  alleles A0/A1 are taken from that id.
* MetaBrain: the table behind the public cis-eQTL browser (MOLGENIS API pages
  saved by src/mr/download_inputs.sh, merged to
  metabrain_cis_eqtl_v2_browser.tsv). Cortex-EUR, QTLRank 1, target genes.
  GRCh38 positions are lifted to GRCh37; the allele frequency comes from the
  1000 Genomes EUR panel (non-palindromic SNVs only; palindromic -> NA).

For both, the rsID is checked against the 1000G EUR LD panel (.bim): a
variant whose rsID is in the panel at a different GRCh37 position is dropped,
and one absent from the panel is kept (it is dropped later by clumping, which
logs it). Every count is logged.
"""

from __future__ import annotations

import csv
import math
import subprocess
from collections import Counter
from pathlib import Path

import yaml

from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.common.runlog import start_run
from src.mr.liftover import DEFAULT_CHAIN, lift
from src.mr.xlsx import iter_rows

RAW = REPO_ROOT / "data" / "raw" / "mr"
OUT = REPO_ROOT / "data" / "interim" / "mr"
UKBPPP_XLSX = RAW / "ukbppp_sun2023" / "41586_2023_6592_MOESM3_ESM.xlsx"
METABRAIN_TSV = RAW / "metabrain_deklein2023" / "metabrain_cis_eqtl_v2_browser.tsv"
BFILE = RAW / "ld_1kg_eur" / "g1000_eur"
EXTRA_GENES = ["HMGCR"]
# Sun et al. 2023, Results: discovery cohort n = 34,557 (data/extracted/mr_qtl_sources.csv).
UKBPPP_DISCOVERY_N = 34557
PALINDROMES = ({"A", "T"}, {"C", "G"})


def target_genes() -> list[str]:
    cfg = yaml.safe_load((CONFIG_DIR / "mr_targets.yaml").read_text(encoding="utf-8"))
    return list(cfg["targets"]) + EXTRA_GENES


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_tsv(path: Path) -> tuple[list[str], list[dict]]:
    with path.open(newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return list(r.fieldnames), list(r)


def bim_index(rsids: set[str]) -> dict[str, tuple[str, int, str, str]]:
    out = {}
    with open(BFILE.with_suffix(".bim")) as fh:
        for line in fh:
            c, snp, _, p, a1, a2 = line.split()
            if snp in rsids:
                out[snp] = (c, int(p), a1.upper(), a2.upper())
    return out


def panel_freq(rsids: set[str], tmp: Path) -> dict[str, dict[str, float]]:
    """1000G EUR allele frequencies: rsid -> {allele: frequency}."""
    tmp.parent.mkdir(parents=True, exist_ok=True)
    (tmp.with_suffix(".snps")).write_text("\n".join(sorted(rsids)) + "\n")
    subprocess.run(["plink", "--bfile", str(BFILE), "--extract", str(tmp.with_suffix(".snps")), "--freq",
                    "--out", str(tmp)], check=True, capture_output=True)
    out = {}
    with open(f"{tmp}.frq") as fh:
        next(fh)
        for line in fh:
            _, snp, a1, a2, maf, _ = line.split()
            out[snp] = {a1.upper(): float(maf), a2.upper(): 1 - float(maf)}
    return out


def ukbppp(genes: set[str], run) -> Path:
    rows = iter_rows(UKBPPP_XLSX, "ST9")
    header = None
    for r in rows:
        if r and r[0].startswith("Variant ID"):
            header = r
            break
    if header is None:
        raise RuntimeError("ST9 header row not found")
    ix = {h: i for i, h in enumerate(header)}
    need = ["Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)", "CHROM", "GENPOS (hg38)", "Assay Target", "rsID",
            "A1FREQ (discovery)", "BETA (discovery, wrt. A1)", "SE (discovery)", "log10(p) (discovery)", "cis/trans",
            "Annotated gene consequence"]
    miss = [h for h in need if h not in ix]
    if miss:
        raise RuntimeError(f"ST9 lacks columns {miss}")
    kept, n = [], Counter()
    for r in rows:
        if len(r) <= ix["Assay Target"] or not r[0]:
            continue
        n["ST9 rows read"] += 1
        if r[ix["Assay Target"]] not in genes:
            continue
        n["rows for target genes"] += 1
        vid = r[ix[need[0]]].split(":")
        chrom37_id, pos37_id, a0, a1 = vid[0], vid[1], vid[2], vid[3]
        lp = float(r[ix["log10(p) (discovery)"]])
        kept.append([r[ix["CHROM"]], r[ix["GENPOS (hg38)"]], r[ix["rsID"]], a1, a0, r[ix["A1FREQ (discovery)"]],
                     r[ix["BETA (discovery, wrt. A1)"]], r[ix["SE (discovery)"]], r[ix["log10(p) (discovery)"]],
                     f"{10 ** (-lp):.6g}" if lp < 300 else "0", str(UKBPPP_DISCOVERY_N), r[ix["Assay Target"]],
                     r[ix["cis/trans"]], r[ix["Annotated gene consequence"]], chrom37_id, pos37_id])
    hdr = ["CHROM", "GENPOS_hg38", "rsID", "A1", "A0", "A1FREQ", "BETA", "SE", "log10p", "P", "N", "Assay_Target",
           "cis_trans", "consequence", "chrom_hg37_in_id", "pos_hg37_in_id"]
    grch38 = OUT / "ukbppp_discovery_st9_targets_grch38.tsv"
    write_tsv(grch38, hdr, kept)
    lifted = OUT / "ukbppp_discovery_st9_targets_lifted.tsv"
    counts = lift(grch38, lifted, "CHROM", "GENPOS_hg38", DEFAULT_CHAIN, lifted.with_name(lifted.name + ".liftover_log.tsv"))
    for k, v in counts.items():
        run.info(f"UKB-PPP liftover: {k}: {v}")
    h, recs = read_tsv(lifted)
    bim = bim_index({x["rsID"] for x in recs})
    out = []
    for x in recs:
        same_id = (x["chrom_grch37"], x["pos_grch37"]) == (x["chrom_hg37_in_id"], x["pos_hg37_in_id"])
        n["lifted position equals the hg37 position in the variant id" if same_id else
          "lifted position DIFFERS from the variant id (dropped)"] += 1
        if not same_id:
            run.drop(stage="ukbppp_build_check", item=x["rsID"], reason="lifted GRCh37 position != hg37 position in ST9 id")
            continue
        b = bim.get(x["rsID"])
        if b and (b[0], b[1]) != (x["chrom_grch37"], int(x["pos_grch37"])):
            n["rsID at another position in 1000G EUR (dropped)"] += 1
            run.drop(stage="ukbppp_rsid_check", item=x["rsID"], reason="rsID at a different GRCh37 position in 1000G EUR")
            continue
        n["rsID in 1000G EUR at the same position" if b else "rsID absent from 1000G EUR (kept; clumping drops it)"] += 1
        out.append([x[c] for c in h])
    path = OUT / "ukbppp_discovery_st9_targets_grch37.tsv"
    write_tsv(path, h, out)
    for k, v in n.items():
        run.info(f"UKB-PPP: {k}: {v}")
    per = Counter((x[h.index("Assay_Target")], x[h.index("cis_trans")]) for x in out)
    for (g, ct), v in sorted(per.items()):
        run.info(f"UKB-PPP rows: {g} {ct}: {v}")
    return path


def metabrain(genes: set[str], run) -> Path:
    h, recs = read_tsv(METABRAIN_TSV)
    n = Counter()
    kept = []
    for x in recs:
        n["browser rows read"] += 1
        if x["GeneSymbol"] not in genes:
            continue
        n["rows for target genes (all regions, all ranks)"] += 1
        if x["Brainregion"] != "Cortex-EUR":
            continue
        if x["QTLRank"] != "1":
            n["Cortex-EUR rank >= 2 (conditional effect; not used)"] += 1
            run.drop(stage="metabrain_rank", item=f"{x['GeneSymbol']} {x['SNP']}", reason="QTLRank >= 2: conditional effect")
            continue
        chrom, pos, rsid, _ = x["SNP"].split(":", 3)
        alleles = x["SNPAlleles"].upper().split("/")
        ea = x["SNPEffectAllele"].upper()
        if len(alleles) != 2 or ea not in alleles:
            n["allele fields inconsistent (dropped)"] += 1
            run.drop(stage="metabrain_alleles", item=rsid, reason="effect allele not one of SNPAlleles")
            continue
        oa = alleles[0] if alleles[1] == ea else alleles[1]
        n["Cortex-EUR rank 1 rows kept"] += 1
        kept.append([x["Brainregion"], x["GeneSymbol"], x["Gene"], chrom, pos, rsid, ea, oa, x["MetaBeta"], x["MetaSE"],
                     x["MetaP"], x["MetaPN"], x["QTLRank"]])
    hdr = ["Brainregion", "GeneSymbol", "Gene", "SNPChr_hg38", "SNPPos_hg38", "rsid", "SNPEffectAllele", "other_allele",
           "MetaBeta", "MetaSE", "MetaP", "MetaPN", "QTLRank"]
    grch38 = OUT / "metabrain_cortex_eur_rank1_targets_grch38.tsv"
    write_tsv(grch38, hdr, kept)
    lifted = OUT / "metabrain_cortex_eur_rank1_targets_lifted.tsv"
    counts = lift(grch38, lifted, "SNPChr_hg38", "SNPPos_hg38", DEFAULT_CHAIN, lifted.with_name(lifted.name + ".liftover_log.tsv"))
    for k, v in counts.items():
        run.info(f"MetaBrain liftover: {k}: {v}")
    h, recs = read_tsv(lifted)
    rs = {x["rsid"] for x in recs}
    bim = bim_index(rs)
    frq = panel_freq(rs, OUT / "tmp_metabrain_panel")
    out = []
    for x in recs:
        b = bim.get(x["rsid"])
        if b is None:
            n["rsID absent from 1000G EUR: no allele frequency (eaf NA; dropped by the reader)"] += 1
            eaf = "NA"
        elif (b[0], b[1]) != (x["chrom_grch37"], int(x["pos_grch37"])):
            n["rsID at another GRCh37 position in 1000G EUR (dropped)"] += 1
            run.drop(stage="metabrain_rsid_check", item=x["rsid"], reason="lifted position != 1000G EUR position")
            continue
        elif {b[2], b[3]} != {x["SNPEffectAllele"], x["other_allele"]}:
            n["alleles differ from 1000G EUR (dropped)"] += 1
            run.drop(stage="metabrain_allele_check", item=x["rsid"], reason="alleles differ from 1000G EUR")
            continue
        elif {x["SNPEffectAllele"], x["other_allele"]} in PALINDROMES:
            n["palindromic: eaf NA (dropped by the reader)"] += 1
            eaf = "NA"
        else:
            n["eaf from 1000G EUR (position and alleles match)"] += 1
            eaf = f"{frq[x['rsid']][x['SNPEffectAllele']]:.6g}"
        out.append([x[c] for c in h] + [eaf])
    path = OUT / "metabrain_cortex_eur_rank1_targets_grch37.tsv"
    write_tsv(path, h + ["eaf_1000g_eur"], out)
    for k, v in n.items():
        run.info(f"MetaBrain: {k}: {v}")
    for x in out:
        run.info("MetaBrain Cortex-EUR rank 1: " + " ".join(f"{c}={v}" for c, v in zip(h + ["eaf_1000g_eur"], x)
                                                          if c in ("GeneSymbol", "rsid", "chrom_grch37", "pos_grch37",
                                                                   "MetaBeta", "MetaSE", "MetaP", "eaf_1000g_eur")))
    return path


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    genes = set(target_genes())
    with start_run(module="mr_prepare_qtl_sources", config_path=CONFIG_DIR / "mr_targets.yaml", seed=None,
                   packages=["pyyaml", "pyliftover"]) as run:
        for p in (UKBPPP_XLSX, METABRAIN_TSV, DEFAULT_CHAIN):
            if not p.exists():
                raise FileNotFoundError(f"{p} missing: run bash src/mr/download_inputs.sh")
        run.record_output(ukbppp(genes, run))
        run.record_output(metabrain(genes, run))
        if not math.isfinite(UKBPPP_DISCOVERY_N):
            raise AssertionError
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
