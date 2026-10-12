# testthat sources helper-*.R before the tests. testthat runs from tests/unit/,
# so walk up to the repository root and load the shared R code from there.
lc_test_root <- normalizePath(getwd())
while (!file.exists(file.path(lc_test_root, "PROJECT_BRIEF.md"))) {
  if (dirname(lc_test_root) == lc_test_root) stop("repository root not found")
  lc_test_root <- dirname(lc_test_root)
}
source(file.path(lc_test_root, "src", "common", "provenance.R"))
source(file.path(lc_test_root, "src", "common", "io.R"))
source(file.path(lc_test_root, "src", "mr", "harmonise.R"))
source(file.path(lc_test_root, "src", "mr", "estimators.R"))
for (f in c("regions", "sumstats", "instruments", "coloc", "steiger", "targets")) {
  source(file.path(lc_test_root, "src", "mr", paste0(f, ".R")))
}
