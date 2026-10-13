validate_meta_positive_control <- function(run, settings, root = lc_repo_root()) {
  path <- file.path(root, settings$extraction$file)
  lc_load_extracted(path, "Claude extraction from published papers", run, root)
  x <- ma_load_extraction(lc_analysis_input(path, "Claude extraction", root), settings)
  x <- x[tolower(x$include) == "true" & x$drug_class == settings$positive_control$drug_class, , drop = FALSE]
  x <- ma_to_log_scale(x, settings$extraction$log_symmetry_tolerance)
  ratio <- x[x$effect_measure %in% settings$pooling$combinable_measures, , drop = FALSE]
  pc <- ma_positive_control(ratio, settings)
  ors <- x[x$effect_measure == "OR" & x$comparator_type == "non_user", , drop = FALSE]
  po <- if (nrow(ors) >= settings$pooling$min_studies) ma_pool(ors, settings, "positive_control_OR") else NULL
  list(passed = isTRUE(pc$passed),
       detail = paste0(pc$detail, if (!is.null(po)) sprintf("; OR group: %.2f [%.2f, %.2f], k = %d",
                                                               exp(po$est), exp(po$ci_lower), exp(po$ci_upper), po$k) else ""),
       table = rbind(cbind(target = "meta_positive_control", measure_group = "HR/RR/IRR", pc$pool),
                     if (!is.null(po)) cbind(target = "meta_positive_control", measure_group = "OR", po)))
}
