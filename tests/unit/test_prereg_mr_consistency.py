"""prereg/mr.md must state exactly what config/mr_targets.yaml and
config/mr_methods.yaml do (and what the validation targets hold).

Each case pairs a config value with the phrase in the pre-registration that
fixes it. A config edit that is not mirrored in the document (or the reverse)
fails here, so the document handed over for freezing cannot drift from the
code that will run.
"""
import csv
import re

import pytest
import yaml

from src.common.config import load_config
from src.common.paths import CONFIG_DIR, REPO_ROOT

DOC = re.sub(r"\s+", " ", (REPO_ROOT / "prereg" / "mr.md").read_text(encoding="utf-8"))
T = load_config(CONFIG_DIR / "mr_targets.yaml")
M = load_config(CONFIG_DIR / "mr_methods.yaml")
HMGCR = yaml.safe_load((REPO_ROOT / "tests/validation/targets/mr_ldl_hmgcr_cad.yaml").read_text(encoding="utf-8"))
NULL = yaml.safe_load((REPO_ROOT / "tests/validation/targets/mr_null_gene_desert.yaml").read_text(encoding="utf-8"))

CASES = [
    # (config, key, expected value, phrase that must appear in the prereg)
    (T, "genome_build", "GRCh37", "Genome build GRCh37 everywhere"),
    (T, "outcome.gwas_catalog_id", "GCST90027158", "GWAS Catalog GCST90027158"),
    (T, "outcome.n_cases", 85934, "85,934 cases"),
    (T, "outcome.n_controls", 401577, "401,577 controls"),
    (T, "replication.release", "R13", "FinnGen release R13, endpoint G6_ALZHEIMER"),
    (T, "replication.finngen_endpoint", "G6_ALZHEIMER", "endpoint G6_ALZHEIMER"),
    (T, "replication.n_cases", 16652, "16,652 cases"),
    (T, "replication.n_controls", 483534, "483,534 controls"),
    (T, "instruments.p_max", 5.0e-8, "p < 5×10⁻⁸"),
    (T, "instruments.maf_min", 0.01, "minor-allele frequency ≥ 0.01"),
    (T, "instruments.f_stat_min", 10, "F > 10"),
    (T, "instruments.ld_r2_max", 0.001, "LD clumping at r² < 0.001 within 10 Mb"),
    (T, "instruments.ld_window_kb", 10000, "within 10 Mb"),
    (T, "window_kb_sensitivity", 1000, "Sensitivity: **±1,000 kb**"),
    (T, "positive_control.p_max", 0.05, "p < 0.05 and the expected sign"),
    (T, "coloc.pp4_min", 0.8, "PP4 ≥ 0.8"),
    (T, "coloc.p1", 1.0e-4, "p1 = p2 = 10⁻⁴"),
    (T, "coloc.p2", 1.0e-4, "p1 = p2 = 10⁻⁴"),
    (T, "coloc.p12", 5.0e-6, "p12 = 5×10⁻⁶"),
    (T, "source_selection.hierarchy", ["pqtl_ukbppp", "pqtl_decode", "metabrain", "eqtlgen"],
     "**UKB-PPP, deCODE, MetaBrain, eQTLGen**"),
    (T, "source_selection.blood_only_source", "eqtlgen", "option (a), blood only"),
    (T, "source_selection.orientation", "level_increasing_allele",
     "pQTL instruments to the protein-increasing allele"),
    (T, "exposures.pqtl_decode.file", None, "not used: the open supplement has cis pQTLs among the targets only for LGMN"),
    (T, "gene_deserts.longest_fraction", 0.03, "the 3% longest of the remaining intervals"),
    (T, "gene_deserts.gene_types", ["protein_coding"], "GENCODE release 19 (GRCh37) protein-coding gene records"),
    (T, "gene_deserts.boundary_gap_min_bp", 250000, "assembly gaps longer than 250 kb"),
    (T, "gene_deserts.variant_margin_kb", 100, "at least 100 kb inside a desert from both flanking genes"),
    (T, "gene_deserts.thin_bp", 50000, "thinned to at least 50,000 bp apart"),
    (M, "harmonise.palindrome_maf_max", 0.42, "both minor-allele frequencies are at most 0.42"),
    (M, "crosscheck.abs_tol", 1.0e-10, "absolute 10⁻¹⁰"),
    (M, "crosscheck.rel_tol", 1.0e-8, "relative 10⁻⁸"),
    (M, "crosscheck.wm_se_rel_tol", 0.15, "within 15% relative"),
    (M, "ivw.model", "multiplicative_random_effects", "IVW with multiplicative random effects"),
    (M, "egger.min_variants", 3, "each needs at least 3 instruments"),
    (M, "weighted_median.min_variants", 3, "each needs at least 3 instruments"),
]


def _get(cfg, key):
    try:
        return cfg.get(key)
    except Exception:  # an explicit null is a value here
        node = cfg.as_dict()
        for part in key.split("."):
            node = node[part]
        return node


@pytest.mark.parametrize("cfg, key, value, phrase", CASES, ids=[c[1] for c in CASES])
def test_prereg_states_what_the_config_does(cfg, key, value, phrase):
    assert _get(cfg, key) == value, f"config {key} = {_get(cfg, key)!r}, but the pre-registration fixes {value!r}"
    assert re.sub(r"\s+", " ", phrase) in DOC, f"prereg/mr.md no longer says: {phrase!r}"


def test_window_and_positive_controls_per_target():
    targets = T.get("targets")
    assert {t["window_kb"] for t in targets.values()} == {100}
    assert "Cis window **±100 kb** of the gene body" in DOC
    signs = {"resting_heart_rate": "resting heart rate", "diastolic_blood_pressure": "diastolic blood pressure",
             "systolic_blood_pressure": "systolic blood pressure", "fasting_glucose": "fasting glucose",
             "sleep_duration": "sleep duration"}
    for gene, t in targets.items():
        pc = t["positive_control"]
        if pc is None:
            assert re.search(rf"\| {gene} \|[^|]*\|[^|]*\|[^|]*\| \*\*none defensible", DOC), gene
        else:
            sign = "+" if pc["expected_sign"] == 1 else "−"
            assert re.search(rf"\| {gene} \|[^|]*\|[^|]*\|[^|]*\| {re.escape(signs[pc['trait']] + ' (' + sign + ')')} \|", DOC), gene


def test_validation_targets_are_the_ones_in_the_document():
    p = HMGCR["published"]
    assert p["doi"] == "10.7554/eLife.73873" and "eLife 10:e73873" in DOC
    assert (round(p["estimate"], 2), p["ci_lower"], p["ci_upper"]) == (1.44, 1.24, 1.68)
    assert "OR **1.44** (95% CI **1.24–1.68**)" in DOC
    assert p["method"] == "ivw" and HMGCR["conversion"]["factor"] == 1
    with (REPO_ROOT / "data/extracted/mr_validation_huang2021_hmgcr.csv").open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert sum(r["role"] == "instrument" for r in rows) == 7 and "The 7 instruments" in DOC
    s = NULL["settings"]
    assert (s["n_sets"], s["variants_per_set"], s["alpha"], s["pass_quantile"]) == (1000, 5, 0.05, 0.99)
    assert "1,000 random sets of 5 desert variants" in DOC and "the 99% binomial quantile" in DOC
    assert NULL["inputs"]["desert_variants"] == T.get("gene_deserts.variants_file")
    assert NULL["inputs"]["genome_build"] == T.get("genome_build")


def test_document_carries_no_status_block_no_placeholders_and_states_no_outcome_read():
    head = DOC[:1500].lower()
    assert "status:" not in head and "{{" not in DOC and "}}" not in DOC
    assert "No row of any Alzheimer's or dementia outcome file" in DOC
