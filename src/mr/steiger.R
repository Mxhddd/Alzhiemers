# Steiger directionality test (Hemani et al. 2017, PLoS Genet): instruments
# should explain more variance in the exposure than in the outcome if the
# causal direction is exposure -> outcome.
#
# Variance explained per variant:
#   quantitative trait (standardised):  r^2 = t^2 / (t^2 + n - 2),  t = beta/se
#   case-control trait (log OR):        r^2 = Vg / (Vg + pi^2/3),
#                                       Vg = 2 p (1 - p) logOR^2
# (the latter is the observed logit-scale approximation; it ignores
# ascertainment and disease prevalence. Its bias is not known in sign for this
# outcome, so Steiger results are reported as supportive, never as a gate).
# Instruments are clumped, so per-variant r^2 are summed.
# Test: Fisher z comparison of two independent correlations.

steiger_r2 <- function(beta, se, eaf, type, n) {
  switch(type,
    quantitative = { t <- beta / se; t^2 / (t^2 + n - 2) },
    case_control = { vg <- 2 * eaf * (1 - eaf) * beta^2; vg / (vg + pi^2 / 3) },
    stop(sprintf("unknown trait type '%s'", type)))
}

steiger_test <- function(h, exposure_type, outcome_type, n_exposure, n_outcome) {
  for (v in c("beta_exp", "se_exp", "eaf_exp", "beta_out", "se_out", "eaf_out")) {
    if (!v %in% names(h) || anyNA(h[[v]])) stop(sprintf("steiger_test needs non-missing %s", v))
  }
  r2_exp <- sum(steiger_r2(h$beta_exp, h$se_exp, h$eaf_exp, exposure_type, n_exposure))
  r2_out <- sum(steiger_r2(h$beta_out, h$se_out, h$eaf_out, outcome_type, n_outcome))
  if (r2_exp >= 1 || r2_out >= 1) stop("summed r^2 >= 1: instruments are not independent or inputs are wrong")
  n_e <- stats::median(n_exposure); n_o <- stats::median(n_outcome)
  zstat <- (atanh(sqrt(r2_exp)) - atanh(sqrt(r2_out))) / sqrt(1 / (n_e - 3) + 1 / (n_o - 3))
  list(r2_exposure = r2_exp, r2_outcome = r2_out,
       correct_direction = r2_exp > r2_out,
       pval = 2 * stats::pnorm(-abs(zstat)))
}
