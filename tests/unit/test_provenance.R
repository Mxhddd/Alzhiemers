# Unit tests for src/common/provenance.R. Run with `make test-r`.
library(testthat)
root <- lc_test_root  # set in helper-root.R
cfg <- file.path(root, "config", "provenance_check.yaml")

test_that("header records commit, config hash, lockfiles, seed and package versions", {
  logs <- withr::local_tempdir()
  run <- lc_start_run("r_unit", cfg, seed = 42L, packages = c("digest", "yaml"), logs_dir = logs)
  lc_finish_run(run, "completed")
  h <- lc_read_header(run$log_path)
  expect_identical(h$git$commit, system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE))
  expect_type(h$git$dirty, "logical")
  expect_identical(h$config$path, "config/provenance_check.yaml")
  expect_identical(h$seed, 42L)
  expect_identical(h$language, "R")
  expect_identical(h$packages$digest, as.character(packageVersion("digest")))
  expect_setequal(names(h$lockfiles), c("uv.lock", "renv.lock"))
})

test_that("config hash is identical to the Python implementation", {
  py <- system2("uv", c("--directory", root, "run", "--frozen", "python", "-m", "src.common.provenance", "sha256", cfg),
                stdout = TRUE)
  expect_identical(lc_sha256(cfg), py[[length(py)]])
})

test_that("a missing package or config raises instead of being recorded as unknown", {
  logs <- withr::local_tempdir()
  expect_error(lc_start_run("r_unit", cfg, seed = NULL, packages = "notARealPackage", logs_dir = logs))
  expect_error(lc_start_run("r_unit", file.path(root, "config", "nope.yaml"), seed = NULL,
                            packages = character(), logs_dir = logs))
  expect_error(lc_start_run("r_unit", cfg, packages = character(), logs_dir = logs), "seed is required")
})

test_that("a failing body writes a FAILED footer and the error propagates", {
  logs <- withr::local_tempdir()
  expect_error(
    lc_with_run("r_unit", cfg, seed = NULL, packages = character(), logs_dir = logs,
                body = function(run) stop("boom")),
    "boom"
  )
  log <- list.files(file.path(logs, "r_unit"), full.names = TRUE)
  expect_length(log, 1)
  expect_true(any(grepl("^status: FAILED", readLines(log))))
})

test_that("drops are counted, logged and reported", {
  logs <- withr::local_tempdir()
  report <- file.path(logs, "dropped.md")
  run <- lc_with_run("r_unit", cfg, seed = NULL, packages = character(), logs_dir = logs,
                     body = function(run) {
                       lc_drop(run, "instruments", "rs1", "F < 10")
                       lc_drop(run, "instruments", "rs2", "F < 10")
                       lc_drop(run, "harmonise", "rs3", "palindromic, MAF near 0.5")
                       lc_write_drop_report(run, report)
                     })
  expect_equal(nrow(run$drops), 3)
  expect_true(any(grepl("Total dropped: \\*\\*3\\*\\*", readLines(report))))
  expect_true(any(grepl("^n_dropped: 3", readLines(run$log_path))))
  expect_error(lc_drop(run, "x", "", "reason"))
})
