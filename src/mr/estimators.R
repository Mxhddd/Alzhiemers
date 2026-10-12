# Two-sample MR estimators, implemented in base R from the published formulas.
#
# Why not call TwoSampleMR / MendelianRandomization? They could not be
# installed in the build environment (CRAN/GitHub R repos blocked). These
# implementations follow the same conventions as TwoSampleMR (cited per
# function), are unit-tested against independent computations, and will be
# cross-checked against TwoSampleMR on identical inputs as soon as it can be
# installed (see tests/unit/test_mr_estimators.R, "cross-implementation").
#
# Input: a harmonised data.frame from mr_harmonise() with columns
# beta_exp, se_exp, beta_out, se_out. Output rows share one schema:
#   method, nsnp, b, se, ci_lower, ci_upper, pval, plus method-specific
#   diagnostics (Q, Q_df, Q_pval, egger_intercept, egger_intercept_se,
#   egger_intercept_pval). Methods that are not estimable for the number of
#   variants return a row with status "not estimable: <reason>" and NA numbers,
#   so they appear in the output table rather than silently vanishing.

.mr_row <- function(method, nsnp, b = NA_real_, se = NA_real_, pval = NA_real_, z = NA_real_,
                    status = "estimated", ...) {
  extra <- list(...)
  row <- data.frame(method = method, nsnp = nsnp, b = b, se = se,
                    ci_lower = b - z * se, ci_upper = b + z * se, pval = pval,
                    Q = NA_real_, Q_df = NA_real_, Q_pval = NA_real_,
                    egger_intercept = NA_real_, egger_intercept_se = NA_real_,
                    egger_intercept_pval = NA_real_, status = status,
                    stringsAsFactors = FALSE)
  for (k in names(extra)) row[[k]] <- extra[[k]]
  row
}

.check_harmonised <- function(h) {
  need <- c("beta_exp", "se_exp", "beta_out", "se_out")
  miss <- setdiff(need, names(h))
  if (length(miss)) stop(sprintf("harmonised data missing columns: %s", paste(miss, collapse = ", ")))
  for (k in need) if (anyNA(h[[k]])) stop(sprintf("harmonised data has NA in %s", k))
  if (any(h$beta_exp == 0)) stop("an instrument has beta_exp == 0; its ratio estimate is undefined")
  invisible(TRUE)
}

# Wald ratio for a single variant. SE first order: se_out / |beta_exp|.
mr_wald_ratio <- function(h, z) {
  .check_harmonised(h)
  if (nrow(h) != 1) stop("mr_wald_ratio takes exactly one variant")
  b <- h$beta_out / h$beta_exp
  se <- h$se_out / abs(h$beta_exp)
  .mr_row("wald_ratio", 1L, b, se, 2 * pnorm(-abs(b / se)), z)
}

# Inverse-variance weighted, multiplicative random effects (TwoSampleMR mr_ivw):
# weighted regression of beta_out on beta_exp through the origin, weights
# 1/se_out^2; SE = fixed-effect SE * max(1, residual SE). Cochran's Q reported.
mr_ivw <- function(h, z) {
  .check_harmonised(h)
  n <- nrow(h)
  if (n < 2) return(.mr_row("ivw", n, status = "not estimable: needs >= 2 variants (use Wald ratio)"))
  w <- 1 / h$se_out^2
  b <- sum(w * h$beta_exp * h$beta_out) / sum(w * h$beta_exp^2)
  se_fixed <- sqrt(1 / sum(w * h$beta_exp^2))
  Q <- sum(w * (h$beta_out - b * h$beta_exp)^2)
  sigma <- sqrt(Q / (n - 1))  # residual standard error of the weighted regression
  se <- se_fixed * max(1, sigma)
  .mr_row("ivw", n, b, se, 2 * pnorm(-abs(b / se)), z,
          Q = Q, Q_df = n - 1, Q_pval = pchisq(Q, n - 1, lower.tail = FALSE))
}

# MR-Egger (Bowden 2015; TwoSampleMR mr_egger_regression): orient every
# variant so beta_exp > 0, then weighted regression with intercept, weights
# 1/se_out^2. SEs divided by min(1, residual SE); t distribution, n-2 df.
mr_egger <- function(h, z, min_variants) {
  .check_harmonised(h)
  n <- nrow(h)
  if (n < min_variants) {
    return(.mr_row("mr_egger", n, status = sprintf("not estimable: needs >= %d variants", min_variants)))
  }
  s <- sign(h$beta_exp)
  x <- abs(h$beta_exp)
  y <- h$beta_out * s
  w <- 1 / h$se_out^2
  fit <- stats::lm(y ~ x, weights = w)
  co <- summary(fit)$coefficients
  sigma <- summary(fit)$sigma
  se_b <- co["x", "Std. Error"] / min(1, sigma)
  se_a <- co["(Intercept)", "Std. Error"] / min(1, sigma)
  b <- unname(co["x", "Estimate"]); a <- unname(co["(Intercept)", "Estimate"])
  Q <- sum(w * stats::residuals(fit)^2)
  .mr_row("mr_egger", n, b, se_b, 2 * pt(abs(b / se_b), df = n - 2, lower.tail = FALSE), z,
          Q = Q, Q_df = n - 2, Q_pval = pchisq(Q, n - 2, lower.tail = FALSE),
          egger_intercept = a, egger_intercept_se = se_a,
          egger_intercept_pval = 2 * pt(abs(a / se_a), df = n - 2, lower.tail = FALSE))
}

# Weighted median of ratio estimates (Bowden et al. 2016, Genet Epidemiol 40:304).
# Interpolated: order the ratio estimates, standardised cumulative weights
# minus half of each weight, and interpolate linearly at 0.5.
mr_weighted_median_point <- function(beta_iv, w) {
  o <- order(beta_iv)
  b <- beta_iv[o]
  wt <- w[o] / sum(w[o])
  p <- cumsum(wt) - 0.5 * wt
  if (p[1] >= 0.5) return(b[1])
  below <- max(which(p < 0.5))
  if (below == length(b)) return(b[length(b)])
  b[below] + (b[below + 1] - b[below]) * (0.5 - p[below]) / (p[below + 1] - p[below])
}

mr_weighted_median <- function(h, z, min_variants, n_boot, weights, seed) {
  .check_harmonised(h)
  n <- nrow(h)
  if (n < min_variants) {
    return(.mr_row("weighted_median", n, status = sprintf("not estimable: needs >= %d variants", min_variants)))
  }
  # Orient every variant so beta_exp > 0. Ratio estimates are unchanged, but
  # the bootstrap draws below then do not depend on which allele was called
  # the effect allele, so the SE is exactly invariant to allele coding.
  s <- sign(h$beta_exp)
  h$beta_exp <- h$beta_exp * s
  h$beta_out <- h$beta_out * s
  beta_iv <- h$beta_out / h$beta_exp
  var_iv <- switch(weights,
    first_order = h$se_out^2 / h$beta_exp^2,
    second_order = h$se_out^2 / h$beta_exp^2 + h$beta_out^2 * h$se_exp^2 / h$beta_exp^4,
    stop(sprintf("unknown weighted_median weights '%s'", weights))
  )
  w <- 1 / var_iv
  b <- mr_weighted_median_point(beta_iv, w)
  # Parametric bootstrap (Bowden 2016): resample both betas from their
  # sampling distributions, keep the original weights. Seeded from config;
  # the caller's RNG state is restored afterwards.
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()))
  set.seed(seed)
  boot <- vapply(seq_len(n_boot), function(i) {
    bx <- stats::rnorm(n, h$beta_exp, h$se_exp)
    by <- stats::rnorm(n, h$beta_out, h$se_out)
    mr_weighted_median_point(by / bx, w)
  }, numeric(1))
  se <- stats::sd(boot)
  .mr_row("weighted_median", n, b, se, 2 * pnorm(-abs(b / se)), z)
}

# All estimable methods for one exposure-outcome pair. `settings` is the
# parsed config/mr_methods.yaml.
mr_run <- function(h, settings) {
  z <- stats::qnorm(1 - (1 - settings$ci_level) / 2)
  if (nrow(h) == 0) stop("no harmonised variants: nothing to estimate")
  if (nrow(h) == 1) return(mr_wald_ratio(h, z))
  rbind(
    mr_ivw(h, z),
    mr_egger(h, z, settings$egger$min_variants),
    mr_weighted_median(h, z, settings$weighted_median$min_variants,
                       settings$weighted_median$n_boot, settings$weighted_median$weights,
                       settings$seed)
  )
}

# Per-variant F statistic, approximated as (beta_exp / se_exp)^2 (one
# instrument, one exposure). Reported alongside every estimate.
mr_f_stat <- function(beta_exp, se_exp) (beta_exp / se_exp)^2
