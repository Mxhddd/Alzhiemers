# Session 0 gate, R side: a run that writes a full provenance header. No analysis.
# Usage: Rscript src/common/provenance_check.R config/provenance_check.yaml
options(warn = 2)  # fail loudly: every warning is an error
source("src/common/provenance.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: Rscript src/common/provenance_check.R <config.yaml>")
cfg <- lc_read_yaml(args[[1]])
if (is.null(cfg$seed)) stop("config has no 'seed'")

run <- lc_with_run("provenance_check_r", args[[1]], seed = as.integer(cfg$seed),
                   packages = c("digest", "yaml", "renv"), body = function(run) {
  h <- lc_read_header(run$log_path)
  lc_log(run, sprintf("header commit=%s config_sha256=%s dirty=%s",
                      h$git$commit, h$config$sha256, h$git$dirty))
})
cat(run$log_path, "\n")
