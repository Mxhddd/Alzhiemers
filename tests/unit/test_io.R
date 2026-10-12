# Unit tests for src/common/io.R (mirrors tests/unit/test_io.py).
library(testthat)
root <- lc_test_root

test_that("sim/ and tests/ paths are refused", {
  expect_error(lc_analysis_input(file.path(root, "sim", "x.csv"), "x", root = root), "simulated")
  expect_error(lc_analysis_input(file.path(root, "tests", "unit", "test_io.R"), "x", root = root), "simulated")
})

test_that("paths outside data/ are refused", {
  expect_error(lc_analysis_input(file.path(root, "README.md"), "x", root = root), "outside data/")
})

test_that("a missing hand-supplied file names who provides it and is not created", {
  p <- file.path(root, "data", "raw", "definitely_absent_file.csv")
  expect_error(lc_analysis_input(p, "the researcher", root = root), "the researcher")
  expect_false(file.exists(p))
})
