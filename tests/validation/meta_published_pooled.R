validate_meta_published_pooled <- function(run, settings, root = lc_repo_root()) {
  t <- lc_load_target("meta_published_pooled", root)
  p <- t$published
  x <- lc_load_extracted(file.path(root, t$inputs$studies), t$inputs$provided_by, run, root)
  pub <- x[x$role == "ma_pooled", , drop = FALSE]
  s <- x[x$role == "ma_component", c("study_id", "estimate", "ci_lower", "ci_upper", "ci_level")]
  if (nrow(pub) != 1) stop("need exactly one ma_pooled row")
  if (nrow(s) < 2) stop("need at least two ma_component rows")
  for (col in c("estimate", "ci_lower", "ci_upper", "ci_level")) s[[col]] <- as.numeric(s[[col]])
  est <- as.numeric(pub$estimate); lo <- as.numeric(pub$ci_lower); hi <- as.numeric(pub$ci_upper)
  s <- ma_to_log_scale(s, settings$extraction$log_symmetry_tolerance)
  fit <- metafor::rma(yi = s$yi, sei = s$sei, method = p$method, test = if (isTRUE(p$knapp_hartung)) "knha" else "z",
                      level = 100 * p$ci_level)
  to_log <- function(v) if (p$reported_as == "ratio") log(v) else v
  ours <- as.numeric(fit$b)
  inside <- to_log(lo) <= ours && ours <= to_log(hi)
  list(passed = inside,
       detail = sprintf("%s pooled (k = %d): ours %.3f [%.3f, %.3f] vs published %.2f [%.2f, %.2f] (ratio scale)",
                        p$method, nrow(s), exp(ours), exp(fit$ci.lb), exp(fit$ci.ub), est, lo, hi),
       table = data.frame(target = t$id, k = nrow(s), method = p$method, our_est = ours, our_ci_lower = fit$ci.lb,
                          our_ci_upper = fit$ci.ub, published_est = to_log(est), published_ci_lower = to_log(lo),
                          published_ci_upper = to_log(hi), abs_difference = abs(ours - to_log(est)),
                          within_published_ci = inside, our_I2 = fit$I2 / 100))
}
