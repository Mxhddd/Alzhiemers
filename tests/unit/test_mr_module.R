# Unit tests for Module 1 (src/mr/*.R beyond harmonise/estimators).
# All numbers are hand-built fixtures for checking arithmetic and logic.
library(testthat)

# ---- coloc -------------------------------------------------------------------
test_that("log ABF equals the log ratio of the two normal densities (Wakefield)", {
  b <- 0.12; se <- 0.03; W <- 0.15
  direct <- dnorm(b, 0, sqrt(se^2 + W^2), log = TRUE) - dnorm(b, 0, se, log = TRUE)
  expect_equal(coloc_labf(b, se, W), direct)
})

brute_force_pp <- function(l1, l2, p1, p2, p12) {
  n <- length(l1)
  h3 <- 0
  for (i in 1:n) for (j in 1:n) if (i != j) h3 <- h3 + exp(l1[i] + l2[j])
  h <- c(1, p1 * sum(exp(l1)), p2 * sum(exp(l2)), p1 * p2 * h3, p12 * sum(exp(l1 + l2)))
  h / sum(h)
}

test_that("coloc posteriors match a brute-force sum over all configurations", {
  d1 <- data.frame(SNP = paste0("s", 1:4), beta = c(0.2, 0.05, 0.01, 0.0), se = 0.04, prior_sd = 0.15)
  d2 <- data.frame(SNP = paste0("s", 1:4), beta = c(0.1, 0.02, 0.05, 0.01), se = 0.05, prior_sd = 0.2)
  res <- coloc_abf(d1, d2, 1e-4, 1e-4, 5e-6)
  bf <- brute_force_pp(coloc_labf(d1$beta, d1$se, 0.15), coloc_labf(d2$beta, d2$se, 0.2), 1e-4, 1e-4, 5e-6)
  expect_equal(unname(res$pp), bf, tolerance = 1e-10)
  expect_equal(sum(res$pp), 1)
})

test_that("coloc separates shared, distinct and absent signals", {
  null <- rep(0, 20); se <- rep(0.02, 20)
  mk <- function(b) data.frame(SNP = paste0("s", 1:20), beta = b, se = se, prior_sd = 0.15)
  b1 <- null; b1[5] <- 0.2
  b2 <- null; b2[5] <- 0.2
  b3 <- null; b3[15] <- 0.2
  expect_gt(coloc_abf(mk(b1), mk(b2), 1e-4, 1e-4, 5e-6)$pp[["PP4"]], 0.9)
  expect_gt(coloc_abf(mk(b1), mk(b3), 1e-4, 1e-4, 5e-6)$pp[["PP3"]], 0.9)
  expect_gt(coloc_abf(mk(b1), mk(null), 1e-4, 1e-4, 5e-6)$pp[["PP1"]], 0.9)
  expect_error(coloc_abf(mk(b1), transform(mk(b1), SNP = paste0("x", 1:20)), 1e-4, 1e-4, 5e-6), "no variants shared")
})

# ---- clumping and instruments -------------------------------------------------
x5 <- data.frame(SNP = paste0("rs", 1:6), pos = c(1000, 2000, 50000, 60000, 9e7, 70000),
                 p = c(1e-10, 1e-12, 1e-9, 1e-8, 1e-9, 1e-9), beta = 0.1, se = 0.01, eaf = 0.3)
ld <- diag(5); dimnames(ld) <- list(paste0("rs", 1:5), paste0("rs", 1:5))
ld["rs1", "rs2"] <- ld["rs2", "rs1"] <- 0.5    # linked: rs2 (smaller p) wins
ld["rs3", "rs4"] <- ld["rs4", "rs3"] <- 0.01   # weakly linked but above 0.001
ld["rs2", "rs5"] <- ld["rs5", "rs2"] <- 0.9    # linked but 90 Mb away: outside window
r2_fun <- function(snps) ld[intersect(snps, rownames(ld)), intersect(snps, rownames(ld)), drop = FALSE]

test_that("greedy clumping keeps leads, honours the window, and logs every removal", {
  cl <- mr_clump(x5, r2_fun, r2_max = 0.001, window_kb = 10000)
  expect_setequal(cl$data$SNP, c("rs2", "rs3", "rs5"))
  reason <- setNames(cl$dropped$reason, cl$dropped$SNP)
  expect_match(reason[["rs1"]], "with rs2")
  expect_match(reason[["rs4"]], "with rs3")
  expect_identical(reason[["rs6"]], "absent from LD reference panel")
})

test_that("instrument filters drop by p, MAF and F with reasons", {
  x <- data.frame(SNP = paste0("rs", 1:4), pos = c(1, 2, 3, 4) * 1e6, eaf = c(0.3, 0.005, 0.3, 0.3),
                  beta = c(0.1, 0.1, 0.1, 0.03), se = c(0.01, 0.01, 0.01, 0.01), p = c(1e-20, 1e-20, 1e-3, 1e-9))
  ident <- function(s) { m <- diag(length(s)); dimnames(m) <- list(s, s); m }
  res <- mr_build_instruments(x, list(p_max = 5e-8, maf_min = 0.01, f_stat_min = 10, ld_r2_max = 0.001, ld_window_kb = 10000), ident)
  expect_identical(res$data$SNP, "rs1")
  reason <- setNames(res$dropped$reason, res$dropped$SNP)
  expect_match(reason[["rs2"]], "MAF")
  expect_match(reason[["rs3"]], "^p >")
  expect_match(reason[["rs4"]], "F <= 10")   # F = (0.03/0.01)^2 = 9
  expect_equal(res$data$f_stat, 100)
})

# ---- eQTLGen Z -> beta --------------------------------------------------------
test_that("Z-score conversion matches hand arithmetic and preserves Z", {
  d <- mr_beta_from_z(5, 0.5, 1000)          # sqrt(2 * .25 * 1025) = sqrt(512.5)
  expect_equal(d$se, 1 / sqrt(512.5))
  expect_equal(d$beta, 5 / sqrt(512.5))
  expect_equal(d$beta / d$se, 5)
})

# ---- Steiger ------------------------------------------------------------------
test_that("Steiger r^2 matches hand arithmetic", {
  expect_equal(steiger_r2(1, 0.1, 0.3, "quantitative", 1000), 100 / (100 + 998))
  expect_equal(steiger_r2(0.1, 0.02, 0.3, "case_control", NA), 0.0042 / (0.0042 + pi^2 / 3))
  expect_error(steiger_r2(1, 1, 0.3, "ordinal", 10), "unknown trait type")
})

test_that("Steiger picks the direction with more variance explained", {
  h <- data.frame(beta_exp = 0.3, se_exp = 0.01, eaf_exp = 0.3, beta_out = 0.02, se_out = 0.01, eaf_out = 0.3)
  st <- steiger_test(h, "quantitative", "case_control", 30000, 100000)
  expect_true(st$correct_direction)
  expect_lt(st$pval, 0.05)
})

# ---- gene regions -------------------------------------------------------------
test_that("gene_region uses the cache, records the release, and stops offline on a miss", {
  sandbox <- withr::local_tempdir()
  calls <- 0
  fake <- function(gene) { calls <<- calls + 1
    data.frame(gene = gene, chrom = "16", start = 55689516, end = 55740140, ensembl_gene_id = "ENSG_TEST",
               ensembl_release = "GRCh37 test", retrieved_utc = "2026-01-01T00:00:00.000000Z") }
  r <- mr_gene_region("SLC6A2", 100, root = sandbox, fetch = fake)
  expect_equal(c(r$start, r$end), c(55689516 - 1e5, 55740140 + 1e5))
  expect_identical(r$ensembl_release, "GRCh37 test")
  mr_gene_region("SLC6A2", 50, root = sandbox, fetch = fake)
  expect_equal(calls, 1)                        # second call served from cache
  expect_error(mr_gene_region("ADRB1", 100, root = sandbox, offline = TRUE), "offline")
  expect_error(mr_gene_region("ADRB1", 100, root = sandbox, build = "GRCh38"), "GRCh37")
})

# ---- summary statistics reader -----------------------------------------------
test_that("sumstats reader maps columns, filters region and gene, derives beta from Z, logs drops", {
  sandbox <- withr::local_tempdir()
  dir.create(file.path(sandbox, "data", "raw"), recursive = TRUE)
  f <- file.path(sandbox, "data", "raw", "eqtl.tsv")
  write.table(data.frame(rsid = c("rs1", "rs2", "rs3", "rs4", "rs5"), chr = c("16", "16", "16", "3", "16"),
                         bp = c(100, 200, 300, 100, 400), A1 = c("a", "C", "G", "T", "A"), A2 = c("G", "T", "A", "C", "G"),
                         Z = c(5, -6, NA, 4, 7), N = 1000, AF = c(0.5, 0.3, 0.2, 0.4, NA), Gene = c("G1", "G1", "G1", "G1", "G2")),
              f, sep = "\t", row.names = FALSE, quote = FALSE)
  cols <- list(SNP = "rsid", chrom = "chr", pos = "bp", effect_allele = "A1", other_allele = "A2",
               z = "Z", n = "N", eaf = "AF", gene = "Gene")
  region <- list(chrom = "16", start = 50, end = 350)
  r <- mr_read_sumstats("data/raw/eqtl.tsv", cols, "GRCh37", "test", region, "G1", sandbox)
  expect_setequal(r$data$SNP, c("rs1", "rs2"))
  expect_identical(r$data$effect_allele[r$data$SNP == "rs1"], "A")
  expect_equal(r$data$beta[r$data$SNP == "rs1"], 5 / sqrt(512.5))
  expect_identical(r$dropped$reason, "Z or N missing")              # rs3; rs4 off-region, rs5 other gene
  expect_identical(attr(r$data, "genome_build"), "GRCh37")
  bad <- cols; bad$z <- "Zscore"
  expect_error(mr_read_sumstats("data/raw/eqtl.tsv", bad, "GRCh37", "test", region, "G1", sandbox), "has no column")
  expect_error(mr_read_sumstats(NULL, cols, "GRCh37", "the researcher", region, "G1", sandbox), "the researcher")
})

# ---- targets, positive control, end-to-end target analysis ---------------------
test_that("expected signs follow the pre-registered predictions and the shipped config is valid", {
  expect_identical(mr_expected_sign("protective"), 1L)
  expect_identical(mr_expected_sign("harmful"), -1L)
  expect_true(is.na(mr_expected_sign("exploratory")))
  cfg <- lc_read_yaml(file.path(lc_test_root, "config", "mr_targets.yaml"))
  expect_silent(mr_check_targets(cfg$targets))
  expect_setequal(names(cfg$targets), c("SLC6A2", "ADRB1", "ADRA1A", "ADRA2A", "LGMN", "SLC18A2",
                                        "MAOA", "DBH", "CACNA1D", "HCRTR2"))
  bad <- cfg$targets; bad$ADRB1$predicted_lower_activity <- "beneficial"
  expect_error(mr_check_targets(bad), "ADRB1")
})

test_that("concordance labels are honest about intervals that include the null", {
  r <- function(b, lo, hi) data.frame(b = b, ci_lower = lo, ci_upper = hi)
  expect_identical(mr_concordance(r(0.2, 0.1, 0.3), 1L), "concordant")
  expect_identical(mr_concordance(r(0.2, -0.1, 0.5), 1L), "direction as predicted, CI includes null")
  expect_match(mr_concordance(r(-0.2, -0.3, -0.1), 1L), "^DISCORDANT")
  expect_match(mr_concordance(r(0.2, 0.1, 0.3), NA_integer_), "exploratory")
})

settings <- lc_read_yaml(file.path(lc_test_root, "config", "mr_methods.yaml"))
settings$weighted_median$n_boot <- 50
cfg <- list(instruments = list(p_max = 5e-8, maf_min = 0.01, f_stat_min = 10, ld_r2_max = 0.001, ld_window_kb = 10000),
            positive_control = list(p_max = 0.05), outcome = list(type = "case_control", n_cases = 20000, n_controls = 60000),
            coloc = list(p1 = 1e-4, p2 = 1e-4, p12 = 5e-6, prior_sd_quantitative = 0.15, prior_sd_case_control = 0.2, pp4_min = 0.8))
mk_ss <- function(beta, se, snps = paste0("rs", seq_along(beta))) {
  x <- data.frame(SNP = snps, effect_allele = "A", other_allele = "G", beta = beta, se = se, eaf = 0.3,
                  pos = seq_along(beta) * 1e5, p = 2 * pnorm(-abs(beta / se)), n = 30000)
  attr(x, "genome_build") <- "GRCh37"; x
}
ident <- function(s) { m <- diag(length(s)); dimnames(m) <- list(s, s); m }
exposure <- mk_ss(c(0.3, 0.25, 0.2, 0.01), c(0.01, 0.01, 0.01, 0.01))   # rs4 is not an instrument

make_io <- function(pc_beta, outcome_called = function() NULL) list(
  r2_fun = ident,
  exposure = function(s, region, gene) list(data = exposure, dropped = data.frame(SNP = character(), reason = character())),
  positive_control = function(gene, pc, region) list(data = mk_ss(pc_beta, rep(0.01, 4))),
  outcome = function(which, region) {
    outcome_called()
    if (which != "main") return(NULL)
    list(data = mk_ss(c(0.03, 0.025, 0.02, 0.0), rep(0.01, 4)), dropped = data.frame(SNP = character(), reason = character()))
  })
tcfg <- list(predicted_lower_activity = "protective", window_kb = 100,
             positive_control = list(trait = "resting_heart_rate", expected_sign = 1, type = "quantitative"))
region <- list(chrom = "1", start = 1, end = 1e7)

test_that("a target that passes its positive control is analysed end to end", {
  drops <- list()
  res <- mr_analyse_target("ADRB1", tcfg, "eqtlgen", region, make_io(c(0.1, 0.08, 0.07, 0)), cfg, settings,
                           function(stage, item, reason) drops[[length(drops) + 1]] <<- reason, function(m) NULL)
  main <- res$rows[res$rows$analysis == "main", ]
  expect_setequal(main$method, c("ivw", "mr_egger", "weighted_median"))
  prim <- main[main$is_primary, ]
  expect_equal(prim$b, 0.1, tolerance = 0.02)   # outcome betas are 0.1 x exposure betas
  expect_identical(prim$concordance, "concordant")
  expect_true(any(grepl("^p > 5e-08", unlist(drops))))   # rs4 (z = 1) fails the p filter first
  expect_match(main$susie_coloc[1], "NOT RUN")
})

test_that("a failing positive control excludes the target before the AD outcome is read", {
  touched <- FALSE
  io <- make_io(c(-0.1, -0.08, -0.07, 0), outcome_called = function() touched <<- TRUE)
  drops <- character()
  res <- mr_analyse_target("ADRB1", tcfg, "eqtlgen", region, io, cfg, settings,
                           function(stage, item, reason) drops <<- c(drops, reason), function(m) NULL)
  expect_false(touched)
  expect_match(res$rows$status, "^excluded: positive control")
  expect_true(any(grepl("positive control", drops)))
})

test_that("a target with no positive control is excluded, not analysed", {
  touched <- FALSE
  t2 <- tcfg; t2$positive_control <- NULL
  res <- mr_analyse_target("SLC6A2", t2, "eqtlgen", region, make_io(0, function() touched <<- TRUE), cfg, settings,
                           function(...) NULL, function(m) NULL)
  expect_false(touched)
  expect_match(res$rows$status, "no positive control specified")
})

test_that("the primary source is the first source in the hierarchy that passed its positive control", {
  h <- c("pqtl_ukbppp", "pqtl_decode", "metabrain", "eqtlgen")
  gate <- data.frame(
    gene = c("A", "A", "A", "B", "B", "C", "C"), window_kb = 100,
    exposure_source = c("eqtlgen", "metabrain", "pqtl_ukbppp", "eqtlgen", "metabrain", "eqtlgen", "metabrain"),
    passed = c(TRUE, TRUE, FALSE, TRUE, NA, FALSE, FALSE), stringsAsFactors = FALSE)
  s <- mr_select_source(gate, h)
  expect_identical(s$primary_source[s$gene == "A"], "metabrain")     # pQTL failed, brain before blood
  expect_identical(s$eligible_sources[s$gene == "A"], "metabrain;eqtlgen")
  expect_identical(s$primary_source[s$gene == "B"], "eqtlgen")       # NA (not run) is not eligible
  expect_true(is.na(s$primary_source[s$gene == "C"]))                # every source failed: untestable
  expect_false(s$testable[s$gene == "C"])
  expect_error(mr_select_source(transform(gate, exposure_source = "gtex"), h), "not in the hierarchy")
  cfg <- lc_read_yaml(file.path(lc_test_root, "config", "mr_targets.yaml"))
  expect_setequal(cfg$source_selection$hierarchy, names(cfg$exposures))
  expect_identical(cfg$source_selection$hierarchy[1], "pqtl_ukbppp")
})
