# Validation: LDL cholesterol via HMGCR variants -> coronary artery disease.
# Called by run_validation.R inside a logged run. Returns a one-row comparison;
# passing means our estimate (same method as the paper) lies inside the
# published CI. Everything else in the output is reported, not hidden.

validate_mr_ldl_hmgcr_cad <- function(run, settings, root = lc_repo_root()) {
  target <- lc_load_target("mr_ldl_hmgcr_cad", root)
  inp <- target$inputs
  exposure <- lc_read_sumstats(file.path(root, inp$exposure), inp$provided_by, inp$genome_build, root)
  outcome <- lc_read_sumstats(file.path(root, inp$outcome), inp$provided_by, inp$genome_build, root)

  hz <- mr_harmonise(exposure, outcome, settings$harmonise$palindrome_maf_max)
  for (i in seq_len(nrow(hz$dropped))) {
    lc_drop(run, stage = "harmonise", item = hz$dropped$SNP[i], reason = hz$dropped$reason[i])
  }
  lc_log(run, sprintf("harmonised: %s",
                      paste(names(hz$counts), hz$counts, sep = "=", collapse = ", ")))
  f <- mr_f_stat(hz$data$beta_exp, hz$data$se_exp)
  lc_log(run, sprintf("per-variant F: min %.1f, median %.1f", min(f), stats::median(f)))

  est <- mr_run(hz$data, settings)
  row <- est[est$method == target$published$method, , drop = FALSE]
  if (nrow(row) != 1 || row$status != "estimated") {
    stop(sprintf("method '%s' used by the paper could not be estimated here: %s",
                 target$published$method, paste(est$status, collapse = "; ")))
  }
  cmp <- lc_compare_to_target(row, target, settings$ci_level)
  list(comparison = cmp, estimates = est)
}
