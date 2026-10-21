# Runs the Module 1 validation targets that read NO Alzheimer's data, before
# the pre-registration is frozen: today only mr_ldl_hmgcr_cad. The gene-desert
# null test (mr_null_gene_desert) reads the AD outcome by design, so it is not
# run here; `make validate-mr` runs both after the freeze and records the
# module status. This script uses the same harness, settings and sources as
# tests/validation/run_validation.R, writes
# results/tables/validation_<target>.csv inside a logged run, and does NOT
# record a validation status (a partial run cannot pass the module).
#
# Usage: Rscript src/mr/run_validation_non_ad.R

options(warn = 2)
source("src/common/provenance.R")
source("src/common/io.R")
source("tests/validation/harness.R")
for (f in c("src/mr/harmonise.R", "src/mr/estimators.R", "src/mr/regions.R", "src/mr/sumstats.R")) source(f)

non_ad_targets <- "mr_ldl_hmgcr_cad"
root <- lc_repo_root()
config_path <- file.path(root, "config", "mr_methods.yaml")
settings <- lc_read_yaml(config_path)

lc_with_run("validation_mr_non_ad", config_path, seed = settings$seed,
            packages = c("yaml", "digest"), body = function(run) {
  for (id in non_ad_targets) {
    source(file.path(root, "tests", "validation", paste0(id, ".R")))
    out <- get(paste0("validate_", id))(run, settings, root)
    lc_log(run, sprintf("%s %s - %s", if (isTRUE(out$passed)) "PASS" else "FAIL", id, out$detail),
           if (isTRUE(out$passed)) "INFO" else "ERROR")
    path <- file.path(root, "results", "tables", sprintf("validation_%s.csv", id))
    utils::write.csv(out$table, path, row.names = FALSE)
    lc_record_output(run, path)
    message(sprintf("%s: %s - %s", id, if (isTRUE(out$passed)) "PASS" else "FAIL", out$detail))
  }
})
