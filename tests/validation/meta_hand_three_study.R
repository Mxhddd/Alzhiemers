validate_meta_hand_three_study <- function(run, settings, root = lc_repo_root()) {
  t <- lc_load_target("meta_hand_three_study", root)
  yi <- unlist(t$inputs$yi); vi <- unlist(t$inputs$vi); e <- t$expected
  if (length(yi) != 3 || length(vi) != 3) stop("hand example needs exactly three yi and three vi")
  fe <- metafor::rma(yi = yi, vi = vi, method = "FE")
  dl <- metafor::rma(yi = yi, vi = vi, method = "DL")
  got <- c(fe_estimate = as.numeric(fe$b), fe_se = fe$se, q = fe$QE, tau2_dl = dl$tau2,
           re_dl_estimate = as.numeric(dl$b), re_dl_se = dl$se)
  want <- unlist(e[names(got)])
  diff <- abs(got - want)
  ok <- all(diff <= e$tolerance)
  list(passed = ok,
       detail = paste(sprintf("%s: pipeline %.6f vs hand %.6f", names(got), got, want), collapse = "; "),
       table = data.frame(target = t$id, quantity = names(got), pipeline = got, hand = want,
                          abs_difference = diff, tolerance = e$tolerance, within_tolerance = diff <= e$tolerance))
}
