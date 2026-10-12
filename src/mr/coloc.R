# Colocalisation, single-causal-variant ABF method (Giambartolomei et al. 2014,
# PLoS Genet; the coloc.abf model), implemented in base R because the coloc
# package could not be installed in the build environment.
#
# Approximate Bayes factor per variant (Wakefield 2009), with V = SE^2 and
# W = prior SD^2:  log ABF = 0.5 * ( log(1 - r) + r * z^2 ),  r = W / (V + W).
# Posterior probabilities of the five hypotheses:
#   H0 no association, H1 trait 1 only, H2 trait 2 only,
#   H3 both, different causal variants, H4 both, one shared causal variant.
#
# SuSiE-coloc (multiple causal variants) needs susieR and an LD matrix; it is
# not implemented. When a region has >1 independent instrument, the caller
# reports that SuSiE-coloc was required but not run.

.logsum <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
.logdiff <- function(a, b) a + log1p(-exp(b - a))  # log(exp(a) - exp(b)), a > b

coloc_labf <- function(beta, se, prior_sd) {
  z <- beta / se
  r <- prior_sd^2 / (se^2 + prior_sd^2)
  0.5 * (log(1 - r) + r * z^2)
}

coloc_abf <- function(d1, d2, p1, p2, p12) {
  for (d in list(d1, d2)) {
    miss <- setdiff(c("SNP", "beta", "se", "prior_sd"), names(d))
    if (length(miss)) stop(sprintf("coloc input missing columns: %s", paste(miss, collapse = ", ")))
  }
  common <- intersect(d1$SNP, d2$SNP)
  if (length(common) == 0) stop("coloc: no variants shared by the two datasets")
  d1 <- d1[match(common, d1$SNP), ]; d2 <- d2[match(common, d2$SNP), ]
  l1 <- coloc_labf(d1$beta, d1$se, d1$prior_sd)
  l2 <- coloc_labf(d2$beta, d2$se, d2$prior_sd)
  lsum1 <- .logsum(l1); lsum2 <- .logsum(l2); lsum12 <- .logsum(l1 + l2)
  lH <- c(
    H0 = 0,
    H1 = log(p1) + lsum1,
    H2 = log(p2) + lsum2,
    H3 = log(p1) + log(p2) + .logdiff(lsum1 + lsum2, lsum12),
    H4 = log(p12) + lsum12
  )
  pp <- exp(lH - .logsum(lH))
  names(pp) <- paste0("PP", 0:4)
  list(pp = pp, nsnps = length(common))
}
