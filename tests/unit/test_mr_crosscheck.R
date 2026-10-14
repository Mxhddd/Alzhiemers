# Cross-implementation check: the base-R ports (src/mr/estimators.R,
# src/mr/harmonise.R) against TwoSampleMR on identical, SIMULATED inputs
# generated here with fixed seeds. No empirical data. The same comparison runs
# inside the pipeline on every real estimate (src/mr/crosscheck.R).
library(testthat)

settings <- lc_read_yaml(file.path(lc_test_root, "config", "mr_methods.yaml"))

sim_h <- function(n, seed, b = 0.3, pleio_sd = 0, mixed_signs = TRUE) {
  set.seed(seed)
  bx <- stats::runif(n, 0.05, 0.4) * if (mixed_signs) sample(c(-1, 1), n, TRUE) else 1
  sx <- stats::runif(n, 0.005, 0.03)
  so <- stats::runif(n, 0.01, 0.05)
  by <- b * bx + stats::rnorm(n, 0, pleio_sd) + stats::rnorm(n, 0, so)
  data.frame(SNP = paste0("rs", seq_len(n)), beta_exp = bx, se_exp = sx, beta_out = by, se_out = so)
}

cases <- list(list(n = 1, seed = 1), list(n = 2, seed = 2), list(n = 3, seed = 3),
              list(n = 5, seed = 4, pleio_sd = 0.02), list(n = 10, seed = 5),
              list(n = 30, seed = 6, pleio_sd = 0.05), list(n = 12, seed = 7, mixed_signs = FALSE))

test_that("every estimate agrees with TwoSampleMR on simulated data sets of 1 to 30 instruments", {
  skip_if_not_installed("TwoSampleMR")
  for (cs in cases) {
    h <- do.call(sim_h, cs)
    est <- mr_run(h, settings)
    cc <- mr_crosscheck(h, est, settings)
    bad <- cc$report[!cc$report$agree, ]
    expect_true(nrow(bad) == 0,
                info = sprintf("n = %d: %s", cs$n, paste(sprintf("%s %s: %g vs %g", bad$method, bad$quantity,
                                                                   bad$base_r, bad$twosamplemr), collapse = "; ")))
    expect_identical(cc$est$b, est$b)                       # nothing replaced when they agree
  }
})

test_that("point estimates and analytic SEs agree to floating-point precision, not just within tolerance", {
  skip_if_not_installed("TwoSampleMR")
  h <- sim_h(30, 6, pleio_sd = 0.05)
  r <- mr_crosscheck(h, mr_run(h, settings), settings)$report
  exact <- r[!(r$method == "weighted_median" & r$quantity %in% c("se", "pval")) & !is.na(r$abs_difference), ]
  expect_lt(max(exact$abs_difference), 1e-12)
  wm <- r[r$method == "weighted_median" & r$quantity == "se", ]
  expect_lt(abs(wm$base_r / wm$twosamplemr - 1), settings$crosscheck$wm_se_rel_tol)
})

test_that("a disagreement is reported and the TwoSampleMR value is the one used", {
  skip_if_not_installed("TwoSampleMR")
  h <- sim_h(10, 5)
  est <- mr_run(h, settings)
  wrong <- est
  wrong$b[wrong$method == "ivw"] <- wrong$b[wrong$method == "ivw"] + 0.01
  cc <- mr_crosscheck(h, wrong, settings)
  expect_false(cc$report$agree[cc$report$method == "ivw" & cc$report$quantity == "b"])
  expect_equal(cc$est$b[cc$est$method == "ivw"], est$b[est$method == "ivw"], tolerance = 1e-12)
  expect_match(cc$est$status[cc$est$method == "ivw"], "TwoSampleMR value used")
})

# ---- Harmonisation --------------------------------------------------------------
sim_alleles <- function(n, seed) {
  set.seed(seed)
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  pairs <- list(c("A", "G"), c("C", "T"), c("A", "C"), c("G", "T"), c("A", "T"), c("C", "G"))
  ex <- data.frame(SNP = paste0("rs", seq_len(n)), stringsAsFactors = FALSE)
  out <- ex
  for (i in seq_len(n)) {
    p <- pairs[[sample(length(pairs), 1)]]
    f <- if (paste(p, collapse = "") %in% c("AT", "CG") && stats::runif(1) < 0.4) stats::runif(1, 0.42, 0.58)
         else stats::runif(1, 0.05, 0.38)
    ex[i, c("effect_allele", "other_allele", "eaf")] <- list(p[1], p[2], f)
    case <- sample(c("direct", "swapped", "strand", "strand_swapped", "incompatible"), 1, prob = c(4, 4, 2, 2, 1))
    fo <- pmin(pmax(f + stats::rnorm(1, 0, 0.01), 0.01), 0.99)
    a <- switch(case,
      direct = list(p[1], p[2], fo), swapped = list(p[2], p[1], 1 - fo),
      strand = list(comp[[p[1]]], comp[[p[2]]], fo), strand_swapped = list(comp[[p[2]]], comp[[p[1]]], 1 - fo),
      incompatible = { q <- setdiff(c("A", "C", "G", "T"), c(p, comp[[p[1]]], comp[[p[2]]]))
                       if (length(q) < 2) list(p[1], p[2], fo) else list(q[1], q[2], fo) })
    out[i, c("effect_allele", "other_allele", "eaf")] <- a
  }
  ex$beta <- stats::rnorm(n, 0.1, 0.05); ex$se <- 0.01
  out$beta <- stats::rnorm(n, 0.02, 0.01); out$se <- 0.005
  attr(ex, "genome_build") <- "GRCh37"; attr(out, "genome_build") <- "GRCh37"
  list(ex = ex, out = out)
}

test_that("harmonisation keeps the same variants with the same aligned outcome effects as TwoSampleMR (action 2)", {
  skip_if_not_installed("TwoSampleMR")
  d <- sim_alleles(400, 11)
  ours <- mr_harmonise(d$ex, d$out, palindrome_maf_max = settings$harmonise$palindrome_maf_max)
  e <- data.frame(SNP = d$ex$SNP, beta.exposure = d$ex$beta, se.exposure = d$ex$se,
                  effect_allele.exposure = d$ex$effect_allele, other_allele.exposure = d$ex$other_allele,
                  eaf.exposure = d$ex$eaf, exposure = "x", id.exposure = "x", stringsAsFactors = FALSE)
  o <- data.frame(SNP = d$out$SNP, beta.outcome = d$out$beta, se.outcome = d$out$se,
                  effect_allele.outcome = d$out$effect_allele, other_allele.outcome = d$out$other_allele,
                  eaf.outcome = d$out$eaf, outcome = "y", id.outcome = "y", stringsAsFactors = FALSE)
  tsmr <- suppressMessages(TwoSampleMR::harmonise_data(e, o, action = 2))
  kept_t <- tsmr[tsmr$mr_keep, ]
  expect_setequal(ours$data$SNP, kept_t$SNP)
  m <- match(ours$data$SNP, kept_t$SNP)
  expect_equal(ours$data$beta_out, kept_t$beta.outcome[m], tolerance = 1e-12)
  expect_equal(ours$data$beta_exp, kept_t$beta.exposure[m], tolerance = 1e-12)
})
