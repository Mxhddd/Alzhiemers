validate_meta_published_pooled <- function(run, settings, root = lc_repo_root()) {
  t <- lc_load_target("meta_published_pooled", root)
  p <- t$published
  s <- utils::read.csv(lc_analysis_input(file.path(root, t$inputs$studies), t$inputs$provided_by, root),
                       stringsAsFactors = FALSE)
  need <- c("study_id", "estimate", "ci_lower", "ci_upper", "ci_level")
  if (!identical(names(s), need)) stop(sprintf("studies file needs columns %s", paste(need, collapse = ", ")))
  if (anyNA(s)) stop("studies file has missing values")
  s <- ma_to_log_scale(s, settings$extraction$log_symmetry_tolerance)
  fit <- metafor::rma(yi = s$yi, sei = s$sei, method = p$method, test = if (isTRUE(p$knapp_hartung)) "knha" else "z",
                      level = 100 * p$ci_level)
  to_log <- function(v) if (p$reported_as == "ratio") log(v) else v
  ours <- as.numeric(fit$b)
  inside <- to_log(p$ci_lower) <= ours && ours <= to_log(p$ci_upper)
  list(passed = inside,
       detail = sprintf("%s pooled (k = %d): ours %.4f [%.4f, %.4f] vs published %.4f [%.4f, %.4f] (log scale)",
                        p$method, nrow(s), ours, fit$ci.lb, fit$ci.ub, to_log(p$estimate), to_log(p$ci_lower), to_log(p$ci_upper)),
       table = data.frame(target = t$id, k = nrow(s), method = p$method, our_est = ours, our_ci_lower = fit$ci.lb,
                          our_ci_upper = fit$ci.ub, published_est = to_log(p$estimate),
                          published_ci_lower = to_log(p$ci_lower), published_ci_upper = to_log(p$ci_upper),
                          abs_difference = abs(ours - to_log(p$estimate)), within_published_ci = inside))
}
