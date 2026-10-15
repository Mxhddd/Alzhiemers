# Informal real-data check of the Module 1 machinery (NOT the validation
# target): LDL cholesterol instrumented by cis variants at HMGCR -> coronary
# artery disease, run with exactly the instrument rules and estimators of
# Module 1. Statin target biology says the estimate must be POSITIVE (higher
# LDL -> higher CAD risk). This checks that reading, region lookup, clumping,
# harmonisation and estimation work end to end on real summary statistics.
#
# It does NOT replace tests/validation/targets/mr_ldl_hmgcr_cad.yaml: that
# target compares with a published CI that the researcher transcribes, and its
# inputs (data/raw/validation/) are researcher-supplied. No published number is
# used or compared here.
#
# Usage: Rscript src/mr/pipeline_check_hmgcr.R
# Output: results/tables/mr_pipeline_check_hmgcr_ldl_cad.csv

options(warn = 2)
source("src/common/provenance.R")
source("src/common/io.R")
for (f in c("harmonise", "estimators", "crosscheck", "regions", "sumstats", "instruments")) {
  source(file.path("src", "mr", paste0(f, ".R")))
}
root <- lc_repo_root()
cfg_path <- file.path(root, "config", "mr_targets.yaml")
cfg <- lc_read_yaml(cfg_path)
mr_settings <- lc_read_yaml(file.path(root, "config", "mr_methods.yaml"))

# Inputs: cis-region subsets written by src/mr/prepare_inputs.py from the raw
# GWAS Catalog downloads (data/raw/mr/MANIFEST.tsv). Column maps from the files'
# own headers.
ldl <- list(
  file = "data/interim/mr/positive_controls/ldl_cholesterol_cisregions.tsv",  # GLGC 2013 joint LDL (GCST002222), GRCh37 via SNP_hg19
  columns = list(SNP = "rsid", chrom = "chrom_grch37", pos = "pos_grch37", effect_allele = "A1",
                 other_allele = "A2", eaf = "Freq.A1.1000G.EUR", beta = "beta", se = "se", p = "P-value"))
cad <- list(
  file = "data/interim/mr/positive_controls/coronary_artery_disease_cisregions.tsv",  # CARDIoGRAMplusC4D 1000G (GCST003116), hg19
  columns = list(SNP = "markername", chrom = "chrom_grch37", pos = "pos_grch37", effect_allele = "effect_allele",
                 other_allele = "noneffect_allele", eaf = "effect_allele_freq", beta = "beta", se = "se_dcl", p = "p_dcl"))
window_kb <- 100  # the Module 1 primary cis window

lc_with_run("mr_pipeline_check_hmgcr", cfg_path, seed = cfg$seed,
            packages = c("data.table", "biomaRt", "yaml", "digest", "TwoSampleMR"), body = function(run) {
  lc_mark_unverified(run, "informal machinery check; not a validation against a published estimate")
  region <- mr_gene_region("HMGCR", window_kb, run = run, root = root)
  e <- mr_read_sumstats(ldl$file, ldl$columns, cfg$genome_build, "Module 1 data step (GLGC LDL)", region, NULL, root)
  o <- mr_read_sumstats(cad$file, cad$columns, cfg$genome_build, "Module 1 data step (CAD)", region, NULL, root)
  for (i in seq_len(nrow(e$dropped))) lc_drop(run, "exposure_read", e$dropped$SNP[i], e$dropped$reason[i])
  for (i in seq_len(nrow(o$dropped))) lc_drop(run, "outcome_read", o$dropped$SNP[i], o$dropped$reason[i])
  inst <- mr_build_instruments(e$data, cfg$instruments, mr_plink_r2_fun(file.path(root, cfg$instruments$ld_reference_bfile)))
  for (i in seq_len(nrow(inst$dropped))) lc_drop(run, "instruments", inst$dropped$SNP[i], inst$dropped$reason[i])
  lc_log(run, sprintf("HMGCR LDL instruments: %s", paste(inst$data$SNP, collapse = ", ")))
  hz <- mr_harmonise(inst$data, o$data, mr_settings$harmonise$palindrome_maf_max)
  hz$data <- hz$data[hz$data$SNP %in% inst$data$SNP, , drop = FALSE]
  lc_log(run, sprintf("harmonised instruments: %d (%s)", nrow(hz$data), paste(hz$data$action, collapse = ",")))
  est <- mr_run_checked(hz$data, mr_settings)
  out <- data.frame(check = "HMGCR LDL -> CAD (informal)", window_kb = window_kb,
                    instruments = paste(hz$data$SNP, collapse = ";"),
                    f_min = min(mr_f_stat(hz$data$beta_exp, hz$data$se_exp)),
                    est[, c("method", "nsnp", "b", "se", "ci_lower", "ci_upper", "pval", "Q", "Q_pval", "status")],
                    or_per_sd_ldl = exp(est$b), stringsAsFactors = FALSE)
  path <- file.path(root, "results", "tables", "mr_pipeline_check_hmgcr_ldl_cad.csv")
  utils::write.csv(out, path, row.names = FALSE)
  lc_record_output(run, path)
  cc <- attr(est, "crosscheck")
  lc_log(run, sprintf("TwoSampleMR cross-check: %d of %d quantities agree", sum(cc$agree), nrow(cc)))
})
