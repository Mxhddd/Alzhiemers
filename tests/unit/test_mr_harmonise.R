# Harmonisation unit test with hand-built strand flips and palindromes
# (PROJECT_BRIEF.md, Module 1 validation). Every row below is constructed by
# hand; the expected action and outcome-beta sign are worked out in the comment.
library(testthat)

ss <- function(df, build = "GRCh37") { attr(df, "genome_build") <- build; df }

exposure <- ss(data.frame(
  SNP           = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs10", "rs12"),
  effect_allele = c("A",   "A",   "A",   "A",   "A",   "A",   "C",   "A",   "A",   "A",    "A"),
  other_allele  = c("G",   "G",   "G",   "G",   "T",   "T",   "G",   "T",   "G",   "G",    "T"),
  beta          = rep(0.1, 11),
  se            = rep(0.01, 11),
  eaf           = c(0.3,   0.3,   0.3,   0.3,   0.2,   0.2,   0.45,  NA,    0.3,   0.3,    0.2),
  stringsAsFactors = FALSE))

outcome <- ss(data.frame(
  SNP           = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs11", "rs12"),
  effect_allele = c("A",   "G",   "T",   "C",   "A",   "A",   "C",   "A",   "A",   "A",    "T"),
  other_allele  = c("G",   "A",   "C",   "T",   "T",   "T",   "G",   "T",   "C",   "G",    "A"),
  beta          = c(0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,   0.05),
  se            = rep(0.02, 11),
  eaf           = c(0.3,   0.7,   0.3,   0.7,   0.2,   0.8,   0.45,  0.2,   0.3,   0.3,    0.8),
  stringsAsFactors = FALSE))

res <- mr_harmonise(exposure, outcome, palindrome_maf_max = 0.42)
kept <- res$data
act <- setNames(kept$action, kept$SNP)
bout <- setNames(kept$beta_out, kept$SNP)
drop_reason <- setNames(res$dropped$reason, res$dropped$SNP)

test_that("non-palindromic variants are aligned from alleles alone", {
  expect_identical(act[["rs1"]], "direct");              expect_equal(bout[["rs1"]], 0.05)   # A/G vs A/G
  expect_identical(act[["rs2"]], "swapped");             expect_equal(bout[["rs2"]], -0.05)  # A/G vs G/A -> negate
  expect_identical(act[["rs3"]], "strand_flip");         expect_equal(bout[["rs3"]], 0.05)   # A/G vs T/C = A/G on other strand
  expect_identical(act[["rs4"]], "strand_flip_swapped"); expect_equal(bout[["rs4"]], -0.05)  # A/G vs C/T = G/A on other strand
  expect_equal(kept$eaf_out[kept$SNP == "rs2"], 0.3)     # 1 - 0.7
})

test_that("palindromes are resolved by allele frequency only when unambiguous", {
  # rs5: A/T, EAF 0.2 in both -> same orientation.
  expect_identical(act[["rs5"]], "palindrome_direct");  expect_equal(bout[["rs5"]], 0.05)
  # rs6: A/T written identically but outcome EAF 0.8 vs 0.2 -> other strand -> negate.
  expect_identical(act[["rs6"]], "palindrome_flipped"); expect_equal(bout[["rs6"]], -0.05)
  # rs12: outcome written T/A with EAF(T)=0.8 => EAF(A)=0.2, agrees -> just swapped -> negate.
  expect_identical(act[["rs12"]], "palindrome_flipped"); expect_equal(bout[["rs12"]], -0.05)
  # rs7: C/G with MAF 0.45 > 0.42 -> strand cannot be inferred -> dropped.
  expect_match(drop_reason[["rs7"]], "MAF above 0.42")
  # rs8: A/T with exposure EAF missing -> dropped.
  expect_identical(drop_reason[["rs8"]], "palindromic, EAF missing")
})

test_that("incompatible and unmatched variants are dropped with reasons, never silently", {
  expect_identical(drop_reason[["rs9"]], "alleles incompatible")   # A/G vs A/C
  expect_identical(drop_reason[["rs10"]], "absent from outcome")
  expect_identical(drop_reason[["rs11"]], "absent from exposure")
  expect_equal(nrow(kept) + nrow(res$dropped), 12)                  # 12 distinct SNPs in
  expect_equal(unname(res$counts["dropped"]), 5)
})

test_that("the palindrome threshold comes from the caller, and changing it changes the result", {
  loose <- mr_harmonise(exposure, outcome, palindrome_maf_max = 0.5)
  expect_true("rs7" %in% loose$data$SNP)
  expect_error(mr_harmonise(exposure, outcome), "palindrome_maf_max is required")
})

test_that("build mismatch, missing build, duplicates and NAs stop the run", {
  expect_error(mr_harmonise(exposure, ss(outcome, "GRCh38"), 0.42), "genome build mismatch")
  no_build <- outcome; attr(no_build, "genome_build") <- NULL
  expect_error(mr_harmonise(exposure, no_build, 0.42), "genome_build")
  expect_error(mr_harmonise(ss(rbind(exposure, exposure[1, ])), outcome, 0.42), "duplicated")
  bad <- outcome; bad$beta[1] <- NA
  expect_error(mr_harmonise(exposure, bad, 0.42), "missing values in 'beta'")
})
