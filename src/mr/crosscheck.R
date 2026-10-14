# Cross-check of the base-R estimators against TwoSampleMR on identical inputs
# (prereg/mr.md, analysis step 4): every estimate is re-run with TwoSampleMR,
# every difference is reported, and where the two disagree the TwoSampleMR
# value is the one used.
#
# "Agree" means:
#   * b, se, pval, Q and the Egger intercept within crosscheck$abs_tol (absolute)
#     or crosscheck$rel_tol (relative), whichever is looser. The formulas are
#     the same, so only floating-point noise is allowed.
#   * The weighted-median SE is a parametric-bootstrap SD. TwoSampleMR draws
#     its bootstrap samples in a different order from ours, so the two SEs
#     differ by Monte Carlo error. They agree if the relative difference is
#     within crosscheck$wm_se_rel_tol (config/mr_methods.yaml, with the
#     reasoning).
# Methods that are "not estimable" in base R must be NA in TwoSampleMR too.

.tsmr_one <- function(h, method, settings) {
  args <- list(b_exp = h$beta_exp, b_out = h$beta_out, se_exp = h$se_exp, se_out = h$se_out)
  r <- switch(method,
    wald_ratio = do.call(TwoSampleMR::mr_wald_ratio, c(args, list(parameters = list()))),
    ivw = do.call(TwoSampleMR::mr_ivw, c(args, list(parameters = list()))),
    mr_egger = do.call(TwoSampleMR::mr_egger_regression, c(args, list(parameters = list()))),
    weighted_median = do.call(TwoSampleMR::mr_weighted_median,
                              c(args, list(parameters = list(nboot = settings$weighted_median$n_boot)))),
    stop(sprintf("no TwoSampleMR counterpart for method '%s'", method)))
  num <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
  list(b = num(r$b), se = num(r$se), pval = num(r$pval), Q = num(r$Q),
       egger_intercept = num(r$b_i), egger_intercept_se = num(r$se_i), egger_intercept_pval = num(r$pval_i))
}

.close <- function(a, b, abs_tol, rel_tol) {
  if (is.na(a) && is.na(b)) return(TRUE)
  if (is.na(a) || is.na(b)) return(FALSE)
  abs(a - b) <= max(abs_tol, rel_tol * max(abs(a), abs(b)))
}

# h: harmonised data (mr_harmonise()$data); est: mr_run() output for h.
# Returns list(est = estimates actually used, report = one row per method x quantity).
mr_crosscheck <- function(h, est, settings) {
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    stop("TwoSampleMR is required for the pre-registered cross-check (prereg/mr.md step 4) but is not installed")
  }
  cc <- settings$crosscheck
  z <- stats::qnorm(1 - (1 - settings$ci_level) / 2)
  report <- list()
  used <- est
  for (i in seq_len(nrow(est))) {
    m <- est$method[i]
    t <- .tsmr_one(h, m, settings)
    ours <- list(b = est$b[i], se = est$se[i], pval = est$pval[i], Q = est$Q[i],
                 egger_intercept = est$egger_intercept[i], egger_intercept_se = est$egger_intercept_se[i],
                 egger_intercept_pval = est$egger_intercept_pval[i])
    if (m != "mr_egger") ours[c("egger_intercept", "egger_intercept_se", "egger_intercept_pval")] <- NA_real_
    if (m %in% c("wald_ratio", "weighted_median")) { ours$Q <- NA_real_; t$Q <- NA_real_ }
    agree_all <- TRUE
    for (q in names(ours)) {
      rel <- if (m == "weighted_median" && q %in% c("se", "pval")) cc$wm_se_rel_tol else cc$rel_tol
      # the weighted-median p value follows its bootstrap SE, so it gets the same Monte Carlo allowance
      ok <- if (m == "weighted_median" && q == "pval") {
        .close(qnorm(ours[[q]] / 2), qnorm(t[[q]] / 2), cc$abs_tol, rel)
      } else .close(ours[[q]], t[[q]], cc$abs_tol, rel)
      agree_all <- agree_all && ok
      report[[length(report) + 1]] <- data.frame(method = m, nsnp = est$nsnp[i], quantity = q,
        base_r = ours[[q]], twosamplemr = t[[q]], abs_difference = abs(ours[[q]] - t[[q]]),
        agree = ok, stringsAsFactors = FALSE)
    }
    if (!agree_all && !is.na(t$b)) {
      used$b[i] <- t$b; used$se[i] <- t$se; used$pval[i] <- t$pval
      used$ci_lower[i] <- t$b - z * t$se; used$ci_upper[i] <- t$b + z * t$se
      if (!is.na(t$Q)) used$Q[i] <- t$Q
      if (m == "mr_egger") {
        used$egger_intercept[i] <- t$egger_intercept; used$egger_intercept_se[i] <- t$egger_intercept_se
        used$egger_intercept_pval[i] <- t$egger_intercept_pval
      }
      used$status[i] <- "estimated (TwoSampleMR value used: base R disagreed, see cross-check report)"
    }
  }
  list(est = used, report = do.call(rbind, report))
}

# mr_run() followed by the cross-check; the report rides along as an attribute
# so callers that only want estimates keep working.
mr_run_checked <- function(h, settings) {
  est <- mr_run(h, settings)
  cc <- mr_crosscheck(h, est, settings)
  out <- cc$est
  attr(out, "crosscheck") <- cc$report
  out
}
