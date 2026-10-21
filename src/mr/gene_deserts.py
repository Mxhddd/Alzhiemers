"""Gene deserts on GRCh37 and the variant list for the Module 1 null test.

    uv run --frozen python -m src.mr.gene_deserts

Rule (config/mr_targets.yaml `gene_deserts`; definition from Ovcharenko et al.
2005, Genome Res 15:137, see data/extracted/mr_gene_desert_definition_ovcharenko2005.csv):

1. Genes = GENCODE release 19 `gene` records of the configured gene types on
   the autosomes; overlapping or touching gene spans are merged.
2. Each chromosome is cut into segments at its ends and at every gap of a
   boundary type (centromere, telomere, heterochromatin, short arm) or longer
   than boundary_gap_min_bp (NCBI GRCh37.p13 gap table).
3. Intergenic intervals = the stretches between two consecutive merged genes
   that lie in the same segment (no boundary gap between them). Stretches that
   touch a segment boundary (before the first or after the last gene of a
   segment) are excluded, as in the paper. Smaller gaps inside an interval are
   kept and counted in its length.
4. Gene deserts = the longest `longest_fraction` of the intervals (count =
   ceil(fraction x N); ties at the cut-off are kept).
5. Variants: 1000G EUR panel SNVs (single-base A/C/G/T alleles, rsID) with
   MAF >= instruments.maf_min lying at least variant_margin_kb inside a desert
   from both flanking genes; thinned (PLINK --bp-space thin_bp), then pruned
   (PLINK --indep-pairwise <ld_window_kb> kb 1 <ld_r2_max>).

Outputs (git-ignored, rebuilt by this script): the desert table, the variant
list (one column, SNP) and PROVENANCE.yaml with input checksums, the rule and
the counts. Reads no outcome data.
"""

from __future__ import annotations

import csv
import gzip
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml

from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.common.provenance import sha256_file
from src.common.runlog import start_run

AUTOSOMES = [str(i) for i in range(1, 23)]


def read_genes(gtf: Path, gene_types: set[str]) -> dict[str, list[tuple[int, int]]]:
    genes = defaultdict(list)
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if f[2] != "gene":
                continue
            chrom = f[0].removeprefix("chr")
            if chrom not in AUTOSOMES:
                continue
            gt = f[8].split('gene_type "', 1)[1].split('"', 1)[0]
            if gt in gene_types:
                genes[chrom].append((int(f[3]), int(f[4])))
    return genes


def merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[list[int]] = []
    for s, e in sorted(spans):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def read_lengths(report: Path) -> dict[str, int]:
    out = {}
    for line in report.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[1] == "assembled-molecule" and f[0] in AUTOSOMES:
            out[f[0]] = int(f[8])
    return out


def read_boundaries(gaps: Path, report: Path, types: set[str], min_bp: int) -> dict[str, list[tuple[int, int]]]:
    """Boundary gaps per chromosome (1-based inclusive), keyed by chromosome name."""
    acc = {}
    for line in report.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[1] == "assembled-molecule" and f[0] in AUTOSOMES:
            acc[f[4]] = f[0]  # GenBank accession (CM000663.1) -> chromosome
    out = defaultdict(list)
    with gzip.open(gaps, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            a, start, stop, length, gtype = line.rstrip("\n").split("\t")[:5]
            if a in acc and (gtype in types or int(length) > min_bp):
                out[acc[a]].append((int(start), int(stop)))
    return out


def intergenic(genes: list[tuple[int, int]], bounds: list[tuple[int, int]]) -> list[tuple[int, int, int, int]]:
    """(start, end, left_gene_end, right_gene_start) for stretches between two
    consecutive genes with no boundary gap in between."""
    out = []
    for (s1, e1), (s2, e2) in zip(genes, genes[1:]):
        lo, hi = e1 + 1, s2 - 1
        if hi < lo:
            continue
        if any(bs <= hi and be >= lo for bs, be in bounds):
            continue  # a segment boundary lies between the two genes: edge stretches, excluded
        out.append((lo, hi, e1, s2))
    return out


def deserts(intervals: dict[str, list[tuple[int, int, int, int]]], fraction: float) -> tuple[list, int, int]:
    flat = [(hi - lo + 1, c, lo, hi) for c, ivs in intervals.items() for lo, hi, _, _ in ivs]
    n = len(flat)
    k = math.ceil(fraction * n)
    lengths = sorted((x[0] for x in flat), reverse=True)
    cut = lengths[k - 1]
    chosen = sorted((x for x in flat if x[0] >= cut), key=lambda x: (AUTOSOMES.index(x[1]), x[2]))
    return chosen, n, cut


def main() -> int:
    cfg = yaml.safe_load((CONFIG_DIR / "mr_targets.yaml").read_text(encoding="utf-8"))
    g = cfg["gene_deserts"]
    ins = cfg["instruments"]
    root = REPO_ROOT
    gtf, gaps, report = root / g["annotation_file"], root / g["gaps_file"], root / g["assembly_report"]
    bfile = root / ins["ld_reference_bfile"]
    regions_out, variants_out = root / g["regions_file"], root / g["variants_file"]
    regions_out.parent.mkdir(parents=True, exist_ok=True)
    with start_run(module="mr_gene_deserts", config_path=CONFIG_DIR / "mr_targets.yaml", seed=None,
                   packages=["pyyaml"]) as run:
        genes = {c: merge(v) for c, v in read_genes(gtf, set(g["gene_types"])).items()}
        bounds = read_boundaries(gaps, report, set(g["boundary_gap_types"]), int(g["boundary_gap_min_bp"]))
        lengths = read_lengths(report)
        for c in AUTOSOMES:
            # chromosome ends are boundaries too (stretches before the first / after the last gene)
            bounds[c] = bounds[c] + [(1, 1), (lengths[c], lengths[c])]
        run.info(f"genes ({', '.join(g['gene_types'])}) after merging: {sum(len(v) for v in genes.values())}")
        ivs = {c: intergenic(genes[c], bounds[c]) for c in AUTOSOMES}
        chosen, n_int, cut = deserts(ivs, float(g["longest_fraction"]))
        total = sum(x[0] for x in chosen)
        run.info(f"intergenic intervals (both flanks genes, same segment): {n_int}")
        run.info(f"gene deserts (longest {g['longest_fraction']}): {len(chosen)}; shortest {cut} bp; "
                 f"total {total} bp")
        n640 = sum(1 for c in ivs for lo, hi, _, _ in ivs[c] if hi - lo + 1 >= 640_000)
        run.info(f"for comparison, intervals >= 640 kb: {n640}")
        with regions_out.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["chrom", "start", "end", "length_bp", "build"])
            for length, c, lo, hi in chosen:
                w.writerow([c, lo, hi, length, "GRCh37"])
        run.record_output(regions_out)

        margin = int(g["variant_margin_kb"]) * 1000
        tmp = regions_out.parent / "tmp"
        with open(f"{tmp}.ranges", "w") as fh:
            for i, (length, c, lo, hi) in enumerate(chosen):
                if hi - margin >= lo + margin:
                    fh.write(f"{c}\t{lo + margin}\t{hi - margin}\td{i}\n")
        # SNVs with rsIDs inside the cores, MAF filter, then thin and prune.
        snvs = tmp.with_suffix(".snv")
        with open(bfile.with_suffix(".bim")) as fi, open(snvs, "w") as fo:
            for line in fi:
                c, snp, _, p, a1, a2 = line.split()
                if c in AUTOSOMES and snp.startswith("rs") and a1 in "ACGT" and a2 in "ACGT" and len(a1) == len(a2) == 1:
                    fo.write(snp + "\n")
        plink = ["plink", "--bfile", str(bfile), "--extract", "range", f"{tmp}.ranges", "--threads", "1"]
        subprocess.run(plink + ["--maf", str(ins["maf_min"]), "--make-bed", "--out", f"{tmp}_core"],
                       check=True, capture_output=True)
        subprocess.run(["plink", "--bfile", f"{tmp}_core", "--extract", str(snvs), "--bp-space", str(g["thin_bp"]),
                        "--make-bed", "--out", f"{tmp}_thin", "--threads", "1"], check=True, capture_output=True)
        subprocess.run(["plink", "--bfile", f"{tmp}_thin", "--indep-pairwise", f"{ins['ld_window_kb']}kb", "1",
                        str(ins["ld_r2_max"]), "--out", f"{tmp}_prune", "--threads", "1"], check=True, capture_output=True)
        n_core = sum(1 for _ in open(f"{tmp}_core.bim"))
        n_thin = sum(1 for _ in open(f"{tmp}_thin.bim"))
        keep = [x.strip() for x in open(f"{tmp}_prune.prune.in") if x.strip()]
        run.info(f"variants: {n_core} in desert cores with MAF >= {ins['maf_min']}; {n_thin} after "
                 f"thinning to {g['thin_bp']} bp; {len(keep)} after LD pruning (r2 < {ins['ld_r2_max']}, "
                 f"{ins['ld_window_kb']} kb)")
        with variants_out.open("w", newline="") as fh:
            fh.write("SNP\n" + "".join(s + "\n" for s in keep))
        run.record_output(variants_out)
        prov = {
            "built_by": "src/mr/gene_deserts.py",
            "rule": {k: g[k] for k in ("longest_fraction", "annotation", "gene_types", "chromosomes",
                                       "boundary_gap_types", "boundary_gap_min_bp", "variant_margin_kb", "thin_bp")},
            "definition_source": "Ovcharenko I et al. Genome Res 2005;15:137-145, doi 10.1101/gr.3015505",
            "build": "GRCh37",
            "inputs": {str(p.relative_to(root)): sha256_file(p) for p in
                       (gtf, gaps, report, bfile.with_suffix(".bim"))},
            "counts": {"merged_genes": sum(len(v) for v in genes.values()), "intergenic_intervals": n_int,
                       "gene_deserts": len(chosen), "shortest_desert_bp": cut, "total_desert_bp": total,
                       "intervals_ge_640kb": n640, "core_variants_maf": n_core, "thinned": n_thin,
                       "variants_after_pruning": len(keep)},
            "outputs": {str(regions_out.relative_to(root)): sha256_file(regions_out),
                        str(variants_out.relative_to(root)): sha256_file(variants_out)},
        }
        (regions_out.parent / "PROVENANCE.yaml").write_text(yaml.safe_dump(prov, sort_keys=False))
        for f in regions_out.parent.glob("tmp*"):
            f.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
