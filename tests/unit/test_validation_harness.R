# Unit tests for tests/validation/harness.R. The target numbers below are
# made-up arithmetic fixtures for testing the comparison logic only; they are
# not published results and never leave this file.
library(testthat)
source(file.path(lc_test_root, "tests", "validation", "harness.R"))

test_that("a target without its published values is refused, and the shipped LDL-CAD target matches its extraction", {
  # A copy of the shipped target with the citation removed must be refused.
  tmp <- file.path(tempfile(), "tests", "validation", "targets")
  dir.create(tmp, recursive = TRUE)
  t <- lc_read_yaml(file.path(lc_test_root, "tests", "validation", "targets", "mr_ldl_hmgcr_cad.yaml"))
  t$published$citation <- NULL
  yaml::write_yaml(t, file.path(tmp, "mr_ldl_hmgcr_cad.yaml"))
  expect_error(lc_load_target("mr_ldl_hmgcr_cad", dirname(dirname(dirname(tmp)))),
               "incomplete - missing: published.citation")
  # The shipped target (filled 2026-10-21) loads, and its numbers are the
  # extracted ones (data/extracted/mr_validation_huang2021_hmgcr.csv).
  shipped <- lc_load_target("mr_ldl_hmgcr_cad", lc_test_root)
  ex <- utils::read.csv(file.path(lc_test_root, "data", "extracted", "mr_validation_huang2021_hmgcr.csv"),
                        stringsAsFactors = FALSE)
  pub <- ex[ex$role == "published_estimate", ]
  expect_equal(shipped$published$estimate, as.numeric(pub$estimate))
  expect_equal(c(shipped$published$ci_lower, shipped$published$ci_upper), as.numeric(c(pub$ci_lower, pub$ci_upper)))
  expect_identical(shipped$published$method, "ivw")
})

fixture <- list(id = "fixture", conversion = list(factor = 1),
                published = list(reported_as = "ratio", estimate = 2, ci_lower = 1.5, ci_upper = 3, ci_level = 0.95))
row <- data.frame(method = "ivw", nsnp = 4, b = log(2.2), se = 0.1,
                  ci_lower = log(2.2) - 0.196, ci_upper = log(2.2) + 0.196)

test_that("comparison is on the log scale and both directions are reported", {
  cmp <- lc_compare_to_target(row, fixture, 0.95)
  expect_equal(cmp$published_ci_lower_log, log(1.5))
  expect_true(cmp$within_published_ci)          # log 2.2 in [log 1.5, log 3]
  expect_true(cmp$published_within_our_ci)      # log 2 = 0.693 in [0.592, 0.984]
})

test_that("an estimate outside the published CI is reported as a failure", {
  far <- row; far$b <- log(4); far$ci_lower <- far$b - 0.1; far$ci_upper <- far$b + 0.1
  expect_false(lc_compare_to_target(far, fixture, 0.95)$within_published_ci)
})

test_that("the conversion factor is applied, including sign flips", {
  flipped <- fixture; flipped$conversion$factor <- -1
  neg <- row; neg$b <- -log(2.2); neg$ci_lower <- neg$b - 0.196; neg$ci_upper <- neg$b + 0.196
  cmp <- lc_compare_to_target(neg, flipped, 0.95)
  expect_true(cmp$within_published_ci)
  expect_lt(cmp$our_ci_lower_log, cmp$our_ci_upper_log)
})

test_that("end to end on synthetic data in a temporary sandbox: target -> inputs -> harmonise -> MR -> compare", {
  # Everything here is generated in a temp directory and deleted afterwards.
  sandbox <- withr::local_tempdir()
  dir.create(file.path(sandbox, "tests", "validation", "targets"), recursive = TRUE)
  dir.create(file.path(sandbox, "data", "raw", "validation"), recursive = TRUE)
  set.seed(5)
  bx <- c(0.08, 0.12, 0.05, 0.10, 0.07)
  exposure <- data.frame(SNP = paste0("rs", 1:5), effect_allele = c("A", "C", "G", "A", "T"),
                         other_allele = c("G", "T", "A", "C", "C"), beta = bx, se = 0.005,
                         eaf = c(0.3, 0.2, 0.4, 0.1, 0.35), source = "synthetic")
  outcome <- exposure
  outcome$beta <- 0.5 * bx + rnorm(5, 0, 0.002); outcome$se <- 0.004
  outcome[2, c("effect_allele", "other_allele")] <- c("T", "C"); outcome$beta[2] <- -outcome$beta[2]  # swapped
  outcome$eaf[2] <- 0.8
  write.csv(exposure, file.path(sandbox, "data/raw/validation/x_exp.csv"), row.names = FALSE)
  write.csv(outcome, file.path(sandbox, "data/raw/validation/x_out.csv"), row.names = FALSE)
  yaml::write_yaml(list(
    id = "fixture", module = "mr", kind = "published_estimate",
    published = list(citation = "synthetic", doi = "none", locator = "none", method = "ivw",
                     scale = "log", reported_as = "log", estimate = 0.5, ci_lower = 0.45,
                     ci_upper = 0.55, ci_level = 0.95, transcribed_by = "test", transcribed_on = "2026-01-01"),
    conversion = list(factor = 1, reason = "same scale"),
    inputs = list(genome_build = "GRCh37", exposure = "data/raw/validation/x_exp.csv",
                  outcome = "data/raw/validation/x_out.csv", provided_by = "test")),
    file.path(sandbox, "tests/validation/targets/fixture.yaml"))

  t <- lc_load_target("fixture", sandbox)
  e <- lc_read_sumstats(file.path(sandbox, t$inputs$exposure), "test", t$inputs$genome_build, sandbox)
  o <- lc_read_sumstats(file.path(sandbox, t$inputs$outcome), "test", t$inputs$genome_build, sandbox)
  hz <- mr_harmonise(e, o, 0.42)
  expect_identical(hz$data$action[hz$data$SNP == "rs2"], "swapped")
  settings <- lc_read_yaml(file.path(lc_test_root, "config", "mr_methods.yaml"))
  est <- mr_run(hz$data, settings)
  cmp <- lc_compare_to_target(est[est$method == "ivw", ], t, settings$ci_level)
  expect_true(cmp$within_published_ci)
  expect_equal(cmp$our_estimate_log, 0.5, tolerance = 0.05)

  # And the sandbox guard still refuses the same file if it sat under sim/.
  dir.create(file.path(sandbox, "sim"))
  file.copy(file.path(sandbox, "data/raw/validation/x_exp.csv"), file.path(sandbox, "sim/x_exp.csv"))
  expect_error(lc_read_sumstats(file.path(sandbox, "sim/x_exp.csv"), "test", "GRCh37", sandbox), "simulated")
})

test_that("the gene-desert null test passes on null data and fails when the outcome is inflated", {
  source(file.path(lc_test_root, "tests", "validation", "mr_null_gene_desert.R"))
  settings <- lc_read_yaml(file.path(lc_test_root, "config", "mr_methods.yaml"))
  run_in <- function(outcome_beta) {
    sandbox <- withr::local_tempdir(.local_envir = parent.frame())
    for (d in c("tests/validation/targets", "data/raw/validation", "config")) dir.create(file.path(sandbox, d), recursive = TRUE)
    n <- 400
    snps <- paste0("rs", seq_len(n))
    write.csv(data.frame(SNP = snps), file.path(sandbox, "data/raw/validation/deserts.csv"), row.names = FALSE)
    write.csv(data.frame(rsid = snps, chr = "1", bp = seq_len(n) * 1e6, ea = "A", oa = "G",
                         b = outcome_beta(n), s = 0.01, af = 0.3),
              file.path(sandbox, "data/raw/ad.csv"), row.names = FALSE)
    yaml::write_yaml(list(genome_build = "GRCh37", outcome = list(file = "data/raw/ad.csv",
      columns = list(SNP = "rsid", chrom = "chr", pos = "bp", effect_allele = "ea", other_allele = "oa",
                     beta = "b", se = "s", eaf = "af"))), file.path(sandbox, "config/mr_targets.yaml"))
    yaml::write_yaml(list(id = "mr_null_gene_desert", module = "mr", kind = "null_test",
      settings = list(seed = 1, n_sets = 300, variants_per_set = 5, alpha = 0.05, pass_quantile = 0.99),
      inputs = list(genome_build = "GRCh37", desert_variants = "data/raw/validation/deserts.csv", provided_by = "test")),
      file.path(sandbox, "tests/validation/targets/mr_null_gene_desert.yaml"))
    run <- lc_start_run("unit_null", file.path(lc_test_root, "config", "mr_methods.yaml"), seed = NULL,
                        packages = character(), logs_dir = file.path(sandbox, "logs"))
    on.exit(lc_finish_run(run, "completed"))
    validate_mr_null_gene_desert(run, settings, sandbox)
  }
  set.seed(11)
  expect_true(run_in(function(n) rnorm(n, 0, 0.01))$passed)        # true null
  expect_false(run_in(function(n) rnorm(n, 0.03, 0.01))$passed)    # every variant "associated"
})
