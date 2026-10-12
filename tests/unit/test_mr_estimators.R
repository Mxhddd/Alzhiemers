# Unit tests for src/mr/estimators.R.
#
# Three kinds of check, none of which uses empirical data:
#   1. hand arithmetic on tiny examples (numbers written out in comments);
#   2. agreement with an independent route to the same quantity (stats::lm);
#   3. statistical behaviour on SIMULATED data generated right here, with the
#      seed fixed in this file: CI coverage, and robustness of the weighted
#      median to invalid instruments (the property Bowden 2016 proved).
library(testthat)

settings <- lc_read_yaml(file.path(lc_test_root, "config", "mr_methods.yaml"))
z95 <- qnorm(0.975)

h3 <- data.frame(beta_exp = c(0.10, 0.20, 0.30), se_exp = c(0.01, 0.01, 0.01),
                 beta_out = c(0.06, 0.09, 0.16), se_out = c(0.02, 0.03, 0.04))

test_that("Wald ratio matches hand arithmetic", {
  r <- mr_wald_ratio(h3[2, ], z95)
  expect_equal(r$b, 0.09 / 0.20)    # 0.45
  expect_equal(r$se, 0.03 / 0.20)   # 0.15
  expect_equal(r$ci_lower, 0.45 - z95 * 0.15)
})

test_that("IVW matches hand arithmetic and a weighted regression through the origin", {
  w <- 1 / h3$se_out^2                            # 2500, 1111.1, 625
  num <- sum(w * h3$beta_exp * h3$beta_out)       # 15 + 20 + 30 = 65
  den <- sum(w * h3$beta_exp^2)                   # 25 + 44.444 + 56.25 = 125.694
  r <- mr_ivw(h3, z95)
  expect_equal(r$b, num / den)
  expect_equal(num, 65, tolerance = 1e-12)
  fit <- lm(beta_out ~ 0 + beta_exp, data = h3, weights = 1 / se_out^2)
  expect_equal(r$b, unname(coef(fit)[1]))
  sigma <- summary(fit)$sigma
  se_fixed <- summary(fit)$coefficients[1, 2] / sigma
  expect_equal(r$se, se_fixed * max(1, sigma))
  expect_equal(r$Q, sum(w * residuals(fit)^2))
})

test_that("IVW never reports an SE below the fixed-effect SE (under-dispersion)", {
  h <- data.frame(beta_exp = c(0.1, 0.2, 0.3), se_exp = 0.01,
                  beta_out = 0.5 * c(0.1, 0.2, 0.3) + c(1e-4, -1e-4, 5e-5), se_out = 0.05)
  r <- mr_ivw(h, z95)
  expect_equal(r$se, sqrt(1 / sum(h$beta_exp^2 / h$se_out^2)))
})

test_that("MR-Egger matches a weighted regression after orienting beta_exp > 0", {
  h <- h3
  h$beta_exp[2] <- -h$beta_exp[2]; h$beta_out[2] <- -h$beta_out[2]   # recode rs2's alleles
  r <- mr_egger(h, z95, min_variants = 3)
  fit <- lm(beta_out ~ beta_exp, data = h3, weights = 1 / se_out^2)  # h3 is already oriented
  expect_equal(r$b, unname(coef(fit)[2]))
  expect_equal(r$egger_intercept, unname(coef(fit)[1]))
  sigma <- summary(fit)$sigma
  expect_equal(r$se, summary(fit)$coefficients[2, 2] / min(1, sigma))
})

test_that("every estimate is invariant to which allele is called the effect allele", {
  flip <- h3
  flip$beta_exp[c(1, 3)] <- -flip$beta_exp[c(1, 3)]
  flip$beta_out[c(1, 3)] <- -flip$beta_out[c(1, 3)]
  a <- mr_run(h3, settings); b <- mr_run(flip, settings)
  expect_equal(a$b, b$b); expect_equal(a$se, b$se)
})

test_that("weighted median point estimate matches hand arithmetic", {
  # Equal weights on 1,2,3: cumulative-minus-half = 1/6, 1/2, 5/6 -> exactly 2.
  expect_equal(mr_weighted_median_point(c(3, 1, 2), c(1, 1, 1)), 2)
  # Equal weights on 1..4: p = .125 .375 .625 .875 -> 2 + (0.5-.375)/.25 = 2.5.
  expect_equal(mr_weighted_median_point(1:4, rep(1, 4)), 2.5)
  # Weights 3,1,1,1: p = .25, .5833, .75, .9167 -> 1 + (0.5-.25)/(.5833-.25) = 1.75.
  expect_equal(mr_weighted_median_point(1:4, c(3, 1, 1, 1)), 1.75)
  # Weights 10,1,1 (=10/12,1/12,1/12): p = 5/12, 7/8, 23/24. The heavy variant
  # does NOT become the median under Bowden's interpolation:
  # 1 + (1/2 - 5/12) / (7/8 - 5/12) = 1 + (1/12)/(11/24) = 1 + 2/11.
  expect_equal(mr_weighted_median_point(1:3, c(10, 1, 1)), 1 + 2 / 11)
  # Only a single variant (cumulative-minus-half = 0.5) returns itself.
  expect_equal(mr_weighted_median_point(7, 1), 7)
})

test_that("weighted median bootstrap is reproducible and leaves the caller's RNG alone", {
  set.seed(1); before <- runif(1); set.seed(1)
  a <- mr_weighted_median(h3, z95, 3, 200, "second_order", seed = 99)
  after <- runif(1)
  expect_identical(before, after)
  b <- mr_weighted_median(h3, z95, 3, 200, "second_order", seed = 99)
  expect_identical(a$se, b$se)
})

test_that("too few variants gives a 'not estimable' row, not a number", {
  r <- mr_egger(h3[1:2, ], z95, min_variants = 3)
  expect_match(r$status, "not estimable")
  expect_true(is.na(r$b))
  expect_equal(nrow(mr_run(h3[1:2, ], settings)), 3)  # all three methods still listed
  expect_identical(mr_run(h3[1, ], settings)$method, "wald_ratio")
})

test_that("bad input stops the run", {
  bad <- h3; bad$se_out[1] <- NA
  expect_error(mr_ivw(bad, z95), "NA in se_out")
  zero <- h3; zero$beta_exp[1] <- 0
  expect_error(mr_ivw(zero, z95), "beta_exp == 0")
})

# ---- Simulated data (generated here; seed fixed in this file) ----------------
sim_instruments <- function(n_snp, b_true, pleio = rep(0, n_snp)) {
  bx_true <- runif(n_snp, 0.05, 0.2)
  se_x <- rep(0.01, n_snp); se_y <- rep(0.02, n_snp)
  data.frame(beta_exp = rnorm(n_snp, bx_true, se_x), se_exp = se_x,
             beta_out = rnorm(n_snp, b_true * bx_true + pleio, se_y), se_out = se_y)
}

test_that("IVW 95% CI covers the true effect about 95% of the time with valid instruments", {
  set.seed(20261013)
  covered <- replicate(400, {
    r <- mr_ivw(sim_instruments(20, b_true = 0.3), z95)
    r$ci_lower <= 0.3 && 0.3 <= r$ci_upper
  })
  expect_gt(mean(covered), 0.92)
  expect_lt(mean(covered), 0.98)
})

test_that("weighted median resists 30% invalid instruments better than IVW (Bowden 2016)", {
  set.seed(20261014)
  err <- replicate(100, {
    pleio <- c(rep(0.05, 6), rep(0, 14))  # 6 of 20 variants with a direct effect on the outcome
    h <- sim_instruments(20, b_true = 0.3, pleio = pleio)
    c(ivw = abs(mr_ivw(h, z95)$b - 0.3),
      wm = abs(mr_weighted_median(h, z95, 3, 50, "second_order", seed = 1)$b - 0.3))
  })
  expect_lt(median(err["wm", ]), median(err["ivw", ]))
})
