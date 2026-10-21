# Per-target drug-target MR: config checks, the positive-control gate, and the
# analysis of one target against the AD outcome.
#
# Order matters and is deliberate: the positive control is evaluated BEFORE
# the AD outcome is touched, so whether a target is excluded can never depend
# on its Alzheimer's result.

.allowed_predictions <- c("protective", "harmful", "exploratory")

# Expected sign of the MR estimate (per higher expression, on log-OR of AD).
# Lower activity protective  => higher activity harmful => positive estimate.
mr_expected_sign <- function(predicted_lower_activity) {
  switch(predicted_lower_activity, protective = 1L, harmful = -1L, exploratory = NA_integer_,
         stop(sprintf("predicted_lower_activity must be one of %s", paste(.allowed_predictions, collapse = ", "))))
}

mr_check_targets <- function(targets) {
  if (!length(targets)) stop("no targets configured")
  for (g in names(targets)) {
    t <- targets[[g]]
    if (!isTRUE(t$predicted_lower_activity %in% .allowed_predictions)) {
      stop(sprintf("target %s: predicted_lower_activity must be one of %s", g, paste(.allowed_predictions, collapse = ", ")))
    }
    if (!is.numeric(t$window_kb) || t$window_kb < 0) stop(sprintf("target %s: window_kb must be a non-negative number", g))
    pc <- t$positive_control
    if (!is.null(pc)) {
      need <- c("trait", "expected_sign", "type")
      miss <- need[vapply(need, function(k) is.null(pc[[k]]), logical(1))]
      if (length(miss)) stop(sprintf("target %s: positive_control lacks %s", g, paste(miss, collapse = ", ")))
      if (!pc$expected_sign %in% c(-1, 1)) stop(sprintf("target %s: positive_control$expected_sign must be 1 or -1", g))
    }
  }
  invisible(TRUE)
}

# Primary estimate: IVW with >= 2 instruments, Wald ratio with one.
mr_primary_row <- function(est) {
  m <- if ("ivw" %in% est$method && est$status[est$method == "ivw"] == "estimated") "ivw" else "wald_ratio"
  row <- est[est$method == m, , drop = FALSE]
  if (nrow(row) != 1 || row$status != "estimated") stop("no primary MR estimate available")
  row
}

# Positive control: the instrument must move a trait it is known to move, in
# the known direction. Returns passed + reason; never throws for a "fail".
mr_positive_control <- function(instruments, pc_sumstats, pc, mr_settings, p_max) {
  hz <- mr_harmonise(instruments, pc_sumstats, mr_settings$harmonise$palindrome_maf_max)
  if (nrow(hz$data) == 0) {
    return(list(passed = FALSE, reason = "no instrument present in the positive-control GWAS", estimate = NULL, harmonise = hz))
  }
  est <- mr_run_checked(hz$data, mr_settings)
  row <- mr_primary_row(est)
  ok_p <- row$pval < p_max
  ok_sign <- sign(row$b) == pc$expected_sign
  reason <- if (ok_p && ok_sign) "passed" else
    sprintf("failed: %s estimate %.3g (p = %.2g); needs p < %g and sign %+d",
            row$method, row$b, row$pval, p_max, pc$expected_sign)
  list(passed = ok_p && ok_sign, reason = reason, estimate = row, harmonise = hz,
       crosscheck = attr(est, "crosscheck"))
}

# Concordance of the primary estimate with the pre-registered prediction.
mr_concordance <- function(row, expected_sign) {
  if (is.na(expected_sign)) return("not pre-specified (exploratory)")
  excludes_null <- row$ci_lower > 0 || row$ci_upper < 0
  same <- sign(row$b) == expected_sign
  if (same && excludes_null) "concordant"
  else if (same) "direction as predicted, CI includes null"
  else if (excludes_null) "DISCORDANT (opposite to prediction, CI excludes null)"
  else "direction opposite to prediction, CI includes null"
}

# One target x one exposure source. `io` bundles the readers so tests can
# supply synthetic files; `drop(stage, item, reason)` logs exclusions.
mr_analyse_target <- function(gene, tcfg, source_name, region, io, cfg, mr_settings, drop, log) {
  base <- list(gene = gene, exposure_source = source_name,
               predicted_lower_activity = tcfg$predicted_lower_activity)
  excluded <- function(reason) {
    drop("target", sprintf("%s [%s]", gene, source_name), reason)
    list(rows = data.frame(base, method = NA_character_, status = paste("excluded:", reason), stringsAsFactors = FALSE))
  }

  exp_raw <- io$exposure(source_name, region, gene)
  for (i in seq_len(nrow(exp_raw$dropped))) drop("exposure_read", exp_raw$dropped$SNP[i], exp_raw$dropped$reason[i])
  if (nrow(exp_raw$data) == 0) return(excluded("no exposure variants in the cis window"))

  inst <- mr_build_instruments(exp_raw$data, cfg$instruments, io$r2_fun)
  for (i in seq_len(nrow(inst$dropped))) drop("instruments", sprintf("%s %s", gene, inst$dropped$SNP[i]), inst$dropped$reason[i])
  if (nrow(inst$data) == 0) return(excluded("no instrument passed p, MAF, F and clumping filters"))
  log(sprintf("%s [%s]: %d instrument(s), F min %.1f", gene, source_name, nrow(inst$data), min(inst$data$f_stat)))

  # ---- Positive-control gate (before the AD outcome is read) ----
  if (is.null(tcfg$positive_control)) return(excluded("no positive control specified (hard gate; researcher to define)"))
  pc_ss <- io$positive_control(gene, tcfg$positive_control, region)
  pc_res <- mr_positive_control(inst$data, pc_ss$data, tcfg$positive_control, mr_settings, cfg$positive_control$p_max)
  log(sprintf("%s [%s] positive control (%s): %s", gene, source_name, tcfg$positive_control$trait, pc_res$reason))
  if (!pc_res$passed) return(excluded(sprintf("positive control %s %s", tcfg$positive_control$trait, pc_res$reason)))

  # ---- AD outcome ----
  out <- io$outcome("main", region)
  for (i in seq_len(nrow(out$dropped))) drop("outcome_read", out$dropped$SNP[i], out$dropped$reason[i])
  if (identical(region$chrom, "X") && nrow(out$data) == 0) {
    return(excluded("chromosome X absent from the outcome summary statistics"))
  }
  hz <- mr_harmonise(inst$data, out$data, mr_settings$harmonise$palindrome_maf_max)
  for (i in seq_len(nrow(hz$dropped))) drop("harmonise", sprintf("%s %s", gene, hz$dropped$SNP[i]), hz$dropped$reason[i])
  if (nrow(hz$data) == 0) return(excluded("no instrument present in the outcome after harmonisation"))
  est <- mr_run_checked(hz$data, mr_settings)
  primary <- mr_primary_row(est)

  # ---- Colocalisation over the whole cis window ----
  prior_out <- if (cfg$outcome$type == "case_control") cfg$coloc$prior_sd_case_control else cfg$coloc$prior_sd_quantitative
  co <- coloc_abf(data.frame(SNP = exp_raw$data$SNP, beta = exp_raw$data$beta, se = exp_raw$data$se,
                             prior_sd = cfg$coloc$prior_sd_quantitative),
                  data.frame(SNP = out$data$SNP, beta = out$data$beta, se = out$data$se, prior_sd = prior_out),
                  cfg$coloc$p1, cfg$coloc$p2, cfg$coloc$p12)
  susie_note <- if (nrow(inst$data) > 1) "SuSiE-coloc required (>1 independent signal) but NOT RUN: susieR unavailable" else "single signal: ABF coloc sufficient"

  # ---- Steiger ----
  n_exp <- if ("n" %in% names(inst$data)) inst$data$n[match(hz$data$SNP, inst$data$SNP)] else stop("exposure sample size column 'n' is required for Steiger")
  st <- steiger_test(hz$data, "quantitative", cfg$outcome$type, n_exp, cfg$outcome$n_cases + cfg$outcome$n_controls)

  f <- inst$data$f_stat[match(hz$data$SNP, inst$data$SNP)]
  egger <- est[est$method == "mr_egger", , drop = FALSE]
  rows <- data.frame(
    base, est[, c("method", "nsnp", "b", "se", "ci_lower", "ci_upper", "pval", "Q", "Q_pval", "status")],
    f_min = min(f), f_median = stats::median(f),
    egger_intercept = if (nrow(egger)) egger$egger_intercept else NA_real_,
    egger_intercept_pval = if (nrow(egger)) egger$egger_intercept_pval else NA_real_,
    coloc_pp4 = co$pp[["PP4"]], coloc_pp3 = co$pp[["PP3"]], coloc_nsnps = co$nsnps,
    coloc_supported = co$pp[["PP4"]] >= cfg$coloc$pp4_min, susie_coloc = susie_note,
    steiger_correct_direction = st$correct_direction, steiger_pval = st$pval,
    positive_control = sprintf("%s: %s", tcfg$positive_control$trait, pc_res$reason),
    is_primary = est$method == primary$method,
    concordance = ifelse(est$method == primary$method,
                         mr_concordance(primary, mr_expected_sign(tcfg$predicted_lower_activity)), NA_character_),
    analysis = "main", stringsAsFactors = FALSE)
  rownames(rows) <- NULL

  # ---- Sensitivity (proxy-free outcome) and replication (FinnGen) ----
  extra <- list()
  cc_parts <- list(positive_control = data.frame(analysis = "positive_control", pc_res$crosscheck),
                   main = data.frame(analysis = "main", attr(est, "crosscheck")))
  for (which in c("proxy_free", "replication")) {
    ss <- io$outcome(which, region)
    if (is.null(ss)) { log(sprintf("%s [%s] %s: not run (no file supplied)", gene, source_name, which)); next }
    hz2 <- mr_harmonise(inst$data, ss$data, mr_settings$harmonise$palindrome_maf_max)
    for (i in seq_len(nrow(hz2$dropped))) drop(sprintf("harmonise_%s", which), sprintf("%s %s", gene, hz2$dropped$SNP[i]), hz2$dropped$reason[i])
    if (nrow(hz2$data) == 0) { drop(which, gene, "no instrument present after harmonisation"); next }
    e2 <- mr_run_checked(hz2$data, mr_settings)
    cc_parts[[which]] <- data.frame(analysis = which, attr(e2, "crosscheck"))
    p2 <- mr_primary_row(e2)
    extra[[which]] <- data.frame(base, p2[, c("method", "nsnp", "b", "se", "ci_lower", "ci_upper", "pval", "Q", "Q_pval", "status")],
                                 is_primary = TRUE,
                                 concordance = mr_concordance(p2, mr_expected_sign(tcfg$predicted_lower_activity)),
                                 analysis = which, stringsAsFactors = FALSE)
  }
  all_rows <- rows
  for (e in extra) {
    for (col in setdiff(names(all_rows), names(e))) e[[col]] <- NA
    all_rows <- rbind(all_rows, e[, names(all_rows)])
  }
  cc <- do.call(rbind, cc_parts)
  list(rows = all_rows, crosscheck = data.frame(gene = gene, exposure_source = source_name, cc,
                                                stringsAsFactors = FALSE))
}

# Stack per-target results whose columns differ (excluded targets carry few).
mr_bind_rows <- function(parts) {
  cols <- unique(unlist(lapply(parts, names)))
  do.call(rbind, lapply(parts, function(p) { for (c in setdiff(cols, names(p))) p[[c]] <- NA; p[, cols] }))
}

# Primary instrument source per target (prereg/mr.md section 5.2; config
# `source_selection`). `gate` has one row per gene x window x exposure_source
# with a logical `passed` (positive control passed; NA = not run). A source is
# eligible only if it passed; the primary source is the first eligible one in
# `hierarchy`. Sources are never pooled. Decided before any outcome was read.
mr_select_source <- function(gate, hierarchy) {
  need <- c("gene", "window_kb", "exposure_source", "passed")
  miss <- setdiff(need, names(gate))
  if (length(miss)) stop(sprintf("gate table lacks %s", paste(miss, collapse = ", ")))
  unknown <- setdiff(unique(stats::na.omit(gate$exposure_source)), hierarchy)
  if (length(unknown)) stop(sprintf("source(s) not in the hierarchy: %s", paste(unknown, collapse = ", ")))
  keys <- unique(gate[, c("gene", "window_kb")])
  out <- lapply(seq_len(nrow(keys)), function(i) {
    g <- gate[gate$gene == keys$gene[i] & gate$window_kb == keys$window_kb[i] & !is.na(gate$exposure_source), ]
    ok <- g$exposure_source[!is.na(g$passed) & g$passed]
    eligible <- hierarchy[hierarchy %in% ok]
    data.frame(gene = keys$gene[i], window_kb = keys$window_kb[i],
               primary_source = if (length(eligible)) eligible[1] else NA_character_,
               eligible_sources = paste(eligible, collapse = ";"),
               testable = length(eligible) > 0, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
