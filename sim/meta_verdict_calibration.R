# SIMULATED calibration of Module 5's contrast verdict under the null (true log ratio 0).
#   OUTDIR=<dir> Rscript sim/meta_verdict_calibration.R
# For k studies with SEs drawn like the dry-run simulator's (0.08-0.3) and between-study SD tau,
# the share of 4,000 data sets whose REML 95% CI excludes 1 above ("AGAINST" for a below-1
# prediction) and on either side, for three CI methods: z, Knapp-Hartung (pre-registered) and
# the ad hoc Knapp-Hartung variant (never narrower than z). Results: prereg/metaanalysis.md change log.
# Share of data sets whose 95% CI excludes 1 on the "AGAINST" side (above 1) and on either side.
set.seed(20261014)
res <- list()
for (k in c(2, 3, 4, 6, 10)) for (tau in c(0, 0.1)) {
  out <- t(replicate(4000, {
    se <- runif(k, 0.08, 0.3); y <- rnorm(k, 0, sqrt(se^2 + tau^2))
    as.vector(sapply(c("z", "knha", "adhoc"), function(tt) {
      f <- tryCatch(metafor::rma(yi = y, sei = se, method = "REML", test = tt), error = function(e) NULL)
      if (is.null(f)) return(c(NA, NA))
      c(f$ci.lb > 0, f$ci.lb > 0 | f$ci.ub < 0)
    }))
  }))
  # out columns: z_against, z_any, knha_against, knha_any, adhoc_against, adhoc_any
  m <- colMeans(out, na.rm = TRUE)
  res[[length(res) + 1]] <- data.frame(k = k, tau = tau, z_against = m[1], z_either = m[2],
    knha_against = m[3], knha_either = m[4], adhoc_against = m[5], adhoc_either = m[6])
}
r <- do.call(rbind, res); rownames(r) <- NULL
write.csv(r, file.path(Sys.getenv("OUTDIR"), "hksj_calibration.csv"), row.names = FALSE); print(format(r, digits = 3))
