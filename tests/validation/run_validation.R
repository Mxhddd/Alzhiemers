# Runs every validation target for one module inside a logged run, writes
# results/tables/validation_<target>.csv, and records the outcome (with the
# code fingerprint) via `python -m src.common.validation record`.
# Usage: Rscript tests/validation/run_validation.R mr
#
# A failing or erroring target is recorded as FAILED with its message; after
# all targets run, any failure stops with a non-zero exit. Nothing is ever
# turned into a pass.

options(warn = 2)  # fail loudly: every warning is an error
source("src/common/provenance.R")
source("src/common/io.R")
source("src/mr/harmonise.R")
source("src/mr/estimators.R")
source("tests/validation/harness.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript tests/validation/run_validation.R <module>")
module <- args[[1]]
root <- lc_repo_root()
config_path <- file.path(root, "config", sprintf("%s_methods.yaml", module))
settings <- lc_read_yaml(config_path)

target_files <- list.files(file.path(root, "tests", "validation", "targets"), pattern = "\\.yaml$", full.names = TRUE)
targets <- Filter(function(f) identical(lc_read_yaml(f)$module, module), target_files)
ids <- sub("\\.yaml$", "", basename(targets))
if (!length(ids)) stop(sprintf("no validation targets for module '%s': an empty validation cannot pass", module))

record_status <- function(results, log_path) {
  passed <- length(results) == length(ids) && all(vapply(results, function(r) isTRUE(r$passed), logical(1)))
  tmp <- tempfile(fileext = ".yaml")
  yaml::write_yaml(list(passed = passed, targets = unname(results)), tmp)
  out <- system2("uv", c("--directory", shQuote(root), "run", "--frozen", "python", "-m",
                         "src.common.validation", "record", module, "--result-file", tmp,
                         "--log", shQuote(lc_relpath(log_path, root))), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(out, "status"))) stop(paste(c("could not record validation status:", out), collapse = "\n"))
  message(paste(out, collapse = "\n"))
  passed
}

lc_with_run(sprintf("validation_%s", module), config_path, seed = settings$seed,
            packages = c("yaml", "digest"), body = function(run) {
  results <- list()
  for (id in ids) {
    source(file.path(root, "tests", "validation", paste0(id, ".R")))
    fn <- get(paste0("validate_", id))
    lc_log(run, sprintf("validation target: %s", id))
    # An error in one target is recorded as that target's FAILURE (with the
    # message) and the remaining targets still run, so every missing input is
    # reported at once. The run as a whole still fails at the end.
    out <- tryCatch(fn(run, settings, root), error = function(e) {
      list(passed = FALSE, detail = paste("ERROR:", conditionMessage(e)),
           table = data.frame(target = id, error = conditionMessage(e)))
    })
    if (!is.logical(out$passed) || length(out$passed) != 1 || is.na(out$passed)) {
      stop(sprintf("validation %s returned no TRUE/FALSE verdict", id))
    }
    lc_log(run, sprintf("%s %s - %s", if (out$passed) "PASS" else "FAIL", id, out$detail),
           if (out$passed) "INFO" else "ERROR")
    table_path <- file.path(root, "results", "tables", sprintf("validation_%s.csv", id))
    dir.create(dirname(table_path), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(out$table, table_path, row.names = FALSE)
    lc_record_output(run, table_path)
    results[[id]] <- list(id = id, passed = out$passed, detail = out$detail)
  }
  if (!record_status(results, run$log_path)) stop(sprintf("validation of module '%s' FAILED", module))
})
