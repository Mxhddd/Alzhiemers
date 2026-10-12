# Homogeneous simulated meta-analyses (generated here, seeded, never stored).
validate_meta_sim_homogeneous <- function(run, settings, root = lc_repo_root()) {
  t <- lc_load_target("meta_sim_homogeneous", root)
  s <- t$settings
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(s$seed)
  i2 <- vapply(seq_len(s$n_reps), function(r) {
    se <- stats::runif(s$k, s$se_min, s$se_max)
    yi <- stats::rnorm(s$k, s$true_log_effect, se)
    ma_rma(yi, se, settings)$I2 / 100
  }, numeric(1))
  med <- stats::median(i2)
  list(passed = med <= s$max_median_i2,
       detail = sprintf("median I^2 %.3f (limit %.3f); %.0f%% of %d replicates have I^2 = 0; 95th percentile %.3f",
                        med, s$max_median_i2, 100 * mean(i2 == 0), s$n_reps, stats::quantile(i2, 0.95)),
       table = data.frame(target = t$id, n_reps = s$n_reps, k = s$k, median_i2 = med,
                          share_zero = mean(i2 == 0), q95_i2 = unname(stats::quantile(i2, 0.95))))
}
