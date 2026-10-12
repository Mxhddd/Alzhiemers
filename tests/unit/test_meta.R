# Unit tests for src/metaanalysis/meta.R. Study numbers are made-up fixtures
# for arithmetic, never published values.
library(testthat)
cfg <- lc_read_yaml(file.path(lc_test_root, "config", "metaanalysis.yaml"))

row <- function(...) {
  base <- list(study_id = "S1", drug_class = "lipophilic_beta_blocker", drug = "propranolol",
               comparator_type = "active", comparator = "atenolol", contrast = "lipophilic_vs_hydrophilic_beta_blocker",
               outcome = "all-cause dementia", effect_measure = "HR", estimate = "0.80", ci_lower = "0.64",
               ci_upper = "1.00", ci_level = "0.95", covariates = "age;sex", brain_penetrance = "high")
  for (d in paste0("robins_d", 1:7)) base[[d]] <- "moderate"
  modifyList(base, list(...))
}
write_extraction <- function(rows) {
  p <- withr::local_tempfile(fileext = ".csv", .local_envir = parent.frame())
  write.csv(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)), p, row.names = FALSE)
  p
}

test_that("a complete extraction loads; any gap or bad label stops with the line number", {
  ok <- ma_load_extraction(write_extraction(list(row(), row(study_id = "S2"))), cfg)
  expect_equal(nrow(ok), 2)
  expect_error(ma_load_extraction(write_extraction(list(row(covariates = ""))), cfg), "line 2: empty covariates")
  expect_error(ma_load_extraction(write_extraction(list(row(comparator_type = "placebo"))), cfg), "comparator_type")
  expect_error(ma_load_extraction(write_extraction(list(row(drug_class = "ssri"))), cfg), "drug_class not listed")
  expect_error(ma_load_extraction(write_extraction(list(row(estimate = "1.20"))), cfg), "outside its own CI")
  expect_error(ma_load_extraction(write_extraction(list(row(robins_d3 = "high"))), cfg), "robins_d3")
  expect_error(ma_load_extraction(write_extraction(list(row(), row())), cfg), "duplicated")
})

test_that("SE is backed out of the CI on the log scale, and asymmetric CIs are refused", {
  x <- data.frame(study_id = "S", estimate = 0.8, ci_lower = 0.64, ci_upper = 1.0, ci_level = 0.95)
  y <- ma_to_log_scale(x, 0.15)
  expect_equal(y$yi, log(0.8))
  expect_equal(y$sei, (log(1.0) - log(0.64)) / (2 * qnorm(0.975)))
  expect_error(ma_to_log_scale(transform(x, estimate = 0.66), 0.15), "not symmetric")
})

test_that("pooling matches hand-computed fixed-effect and DerSimonian-Laird formulas", {
  yi <- c(-0.2, 0.1, -0.4); sei <- c(0.1, 0.2, 0.15)
  w <- 1 / sei^2
  fe <- sum(w * yi) / sum(w)
  Q <- sum(w * (yi - fe)^2)
  tau2 <- max(0, (Q - 2) / (sum(w) - sum(w^2) / sum(w)))
  wr <- 1 / (sei^2 + tau2)
  dl <- sum(wr * yi) / sum(wr)
  x <- data.frame(study_id = c("a", "b", "c"), yi = yi, sei = sei)
  p <- ma_pool(x, cfg, "fixture")
  expect_equal(p$est_dl, dl, tolerance = 1e-8)
  expect_equal(p$Q, Q, tolerance = 1e-8)
  expect_identical(p$status, "pooled")
  expect_lt(p$pi_lower, p$ci_lower)                       # prediction interval wider than CI
})

test_that("a single study is shown with its own transcribed CI, not pooled", {
  x <- data.frame(study_id = "a", yi = log(0.8), sei = 0.1, ci_lower = 0.64, ci_upper = 1.0)
  p <- ma_pool(x, cfg, "one")
  expect_match(p$status, "not pooled: k = 1")
  expect_equal(c(p$ci_lower, p$ci_upper), log(c(0.64, 1.0)))
})

test_that("the positive control passes only when harmful and excluding 1", {
  mk <- function(est) data.frame(study_id = c("a", "b", "c"), drug_class = "anticholinergic",
                                 yi = log(est), sei = c(0.05, 0.06, 0.07))
  expect_true(ma_positive_control(mk(c(1.3, 1.4, 1.25)), cfg)$passed)
  expect_false(ma_positive_control(mk(c(0.95, 1.05, 1.0)), cfg)$passed)
  expect_false(ma_positive_control(mk(c(0.7, 0.75, 0.8)), cfg)$passed)
  expect_match(ma_positive_control(mk(c(1.3, 1.4, 1.25))[0, ], cfg)$detail, "no 'anticholinergic' rows")
})

test_that("triangulation puts drug use, MR and prediction on one 'lower activity' scale and keeps gaps", {
  pooled <- data.frame(drug_class = c("lipophilic_beta_blocker", "guanfacine"), est = c(log(0.8), log(0.9)),
                       ci_lower = c(log(0.7), log(0.85)), ci_upper = c(log(0.9), log(0.95)), k = c(4, 3))
  mr <- data.frame(gene = "ADRB1", analysis = "main", is_primary = TRUE, b = 0.1, ci_lower = 0.05,
                   ci_upper = 0.15, exposure_source = "eqtlgen")
  targets <- lc_read_yaml(file.path(lc_test_root, "config", "mr_targets.yaml"))$targets
  tri <- ma_triangulate(pooled, mr, targets, cfg)
  b1 <- tri[tri$target == "ADRB1", ]
  expect_equal(b1$obs_est, log(0.8))            # beta-blocker = lower activity: sign kept
  expect_equal(b1$mr_est, -0.1)                 # MR per higher expression -> flipped
  expect_identical(b1$obs_vs_prediction, "agree")
  expect_identical(b1$mr_vs_prediction, "agree")
  a2 <- tri[tri$target == "ADRA2A", ]
  expect_equal(a2$obs_est, -log(0.9))           # agonist: lower activity is the reverse of drug use
  expect_identical(a2$obs_vs_prediction, "agree")   # prediction: lower activity harmful
  expect_identical(a2$mr_vs_prediction, "missing")
  expect_equal(nrow(tri), length(targets))      # every target present, gaps included
  expect_true(is.na(tri$obs_est[tri$target == "LGMN"]))
})
