"""Inputs of the HMGCR validation target (tests/validation/targets/mr_ldl_hmgcr_cad.yaml).

    uv run --frozen python -m src.mr.build_hmgcr_validation_inputs

Researcher decision 2026-10-21: reproduce Huang et al. 2021 (eLife 10:e73873),
whose 7 HMGCR instruments are listed in its Supplementary file 1 Table 3
(transcribed in data/extracted/mr_validation_huang2021_hmgcr.csv, role
"instrument"). For those 7 SNPs this script takes the statistics from the two
open GWAS the paper used, as downloaded by src/mr/download_inputs.sh and
cut to the cis regions by src/mr/prepare_inputs.py:

* exposure: GLGC 2013 jointGwasMc_LDL (GWAS Catalog GCST002222), GRCh37;
* outcome: CARDIoGRAMplusC4D 1000 Genomes additive (GCST003116), GRCh37.

It checks that the GLGC alleles, beta, SE and frequency equal the paper's
Table 3 (they are printed from the same file) and stops on any difference.
Writes data/raw/validation/mr_ldl_hmgcr_cad_{exposure,outcome}.csv with
columns SNP, effect_allele, other_allele, beta, se, eaf, source.
"""

from __future__ import annotations

import csv
from pathlib import Path

from src.common.paths import CONFIG_DIR, REPO_ROOT
from src.common.runlog import start_run

EXTRACTED = REPO_ROOT / "data" / "extracted" / "mr_validation_huang2021_hmgcr.csv"
PC = REPO_ROOT / "data" / "interim" / "mr" / "positive_controls"
OUT = REPO_ROOT / "data" / "raw" / "validation"


def rows(path: Path, key: str) -> dict[str, dict]:
    with path.open(newline="") as fh:
        return {r[key]: r for r in csv.DictReader(fh, delimiter="\t")}


def main() -> int:
    with start_run(module="mr_build_hmgcr_validation_inputs", config_path=CONFIG_DIR / "mr_methods.yaml",
                   seed=None, packages=[]) as run:
        with EXTRACTED.open(newline="") as fh:
            inst = [r for r in csv.DictReader(fh) if r["role"] == "instrument"]
        if len(inst) != 7:
            raise RuntimeError(f"expected the paper's 7 instruments, found {len(inst)}")
        ldl = rows(PC / "ldl_cholesterol_cisregions.tsv", "rsid")
        cad = rows(PC / "coronary_artery_disease_cisregions.tsv", "markername")
        exp, out = [], []
        for r in inst:
            snp = r["item"]
            e, o = ldl.get(snp), cad.get(snp)
            if e is None or o is None:
                raise RuntimeError(f"{snp} missing from the {'LDL' if e is None else 'CAD'} GWAS")
            ea, oa = e["A1"].upper(), e["A2"].upper()
            same = (ea, oa) == (r["effect_allele"], r["other_allele"]) and \
                abs(float(e["beta"]) - float(r["beta"])) < 1e-9 and abs(float(e["se"]) - float(r["se"])) < 1e-9 and \
                abs(float(e["Freq.A1.1000G.EUR"]) - float(r["eaf"])) < 1e-9
            if not same:
                raise RuntimeError(f"{snp}: GLGC file {ea}/{oa} {e['beta']} {e['se']} differs from the paper's Table 3")
            run.info(f"{snp}: GLGC statistics equal Huang 2021 Table 3")
            exp.append([snp, ea, oa, e["beta"], e["se"], e["Freq.A1.1000G.EUR"], "GLGC 2013 jointGwasMc_LDL (GCST002222)"])
            out.append([snp, o["effect_allele"].upper(), o["noneffect_allele"].upper(), o["beta"], o["se_dgc"],
                        o["effect_allele_freq"], "CARDIoGRAMplusC4D 1000G additive (GCST003116)"])
        OUT.mkdir(parents=True, exist_ok=True)
        for name, data in (("exposure", exp), ("outcome", out)):
            path = OUT / f"mr_ldl_hmgcr_cad_{name}.csv"
            with path.open("w", newline="") as fh:
                w = csv.writer(fh, lineterminator="\n")
                w.writerow(["SNP", "effect_allele", "other_allele", "beta", "se", "eaf", "source"])
                w.writerows(data)
            run.record_output(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
