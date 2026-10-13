# Module 5 — meta-analysis and triangulation (metafor).
#
# The extraction (data/extracted/meta_extraction.csv) is built from published
# papers by Claude under the sourcing rules: every row names its source
# (source_id -> data/extracted/SOURCES.md), the exact locator (table/row), the
# quantity, a verbatim snippet, verified (true/false; only a human sets true)
# and machine_rechecked (true/false/discrepancy). load_extraction() validates
# it hard and stops on any gap; nothing is imputed, nothing silently dropped.

.ma_required <- c("study_id", "drug_class", "drug", "comparator_type", "comparator", "contrast",
                  "outcome", "effect_measure", "estimate", "ci_lower", "ci_upper", "ci_level",
                  "covariates", "brain_penetrance", paste0("robins_d", 1:7),
                  "source_id", "locator", "quantity", "verbatim", "verified", "machine_rechecked",
                  "invert_ratio", "include", "exclusion_reason")
# not_assessed: risk of bias not yet rated (stated in the pre-registration;
# the risk-of-bias sensitivity analysis cannot run until ratings exist).
.robins_levels <- c("low", "moderate", "serious", "critical", "no_information", "not_assessed")

# Every random-effects fit goes through here so the estimator, the HKSJ
# choice and the convergence controls are the configured ones.
ma_rma <- function(yi, sei, cfg, method = cfg$pooling$method, knha = isTRUE(cfg$pooling$knapp_hartung), ...) {
  metafor::rma(yi = yi, sei = sei, method = method, test = if (knha) "knha" else "z",
               control = cfg$pooling$rma_control, ...)
}

ma_load_extraction <- function(path, cfg) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE, na.strings = character(0), strip.white = TRUE,
                       colClasses = "character")
  if (!identical(names(x), .ma_required)) {
    stop(sprintf("%s: columns must be exactly, in order: %s\nfound: %s", path,
                 paste(.ma_required, collapse = ", "), paste(names(x), collapse = ", ")))
  }
  problems <- character()
  blank <- which(x == "", arr.ind = TRUE)
  if (nrow(blank)) {
    problems <- c(problems, sprintf("line %d: empty %s (write 'none' / 'not_reported' explicitly)",
                                    blank[, 1] + 1, names(x)[blank[, 2]]))
  }
  if (anyDuplicated(x[, c("study_id", "drug", "comparator", "outcome")])) {
    problems <- c(problems, "duplicated (study_id, drug, comparator, outcome) rows")
  }
  num <- c("estimate", "ci_lower", "ci_upper", "ci_level")
  # decimals as printed (for the rounding allowance of the symmetry check)
  ndec <- function(v) ifelse(grepl(".", v, fixed = TRUE), nchar(sub("^[^.]*\\.", "", v)), 0L)
  x$dec_est <- ndec(x$estimate); x$dec_lo <- ndec(x$ci_lower); x$dec_hi <- ndec(x$ci_upper)
  for (col in num) {
    v <- suppressWarnings(as.numeric(x[[col]]))
    bad <- which(is.na(v) & x[[col]] != "")
    if (length(bad)) problems <- c(problems, sprintf("line %d: %s is not a number", bad + 1, col))
    x[[col]] <- v
  }
  chk <- function(cond, msg) { bad <- which(!cond); if (length(bad)) problems <<- c(problems, sprintf("line %d: %s", bad + 1, msg)) }
  if (!length(problems)) {
    chk(x$estimate > 0 & x$ci_lower > 0 & x$ci_upper > 0, "ratio estimates and CI bounds must be > 0")
    chk(x$ci_lower <= x$estimate & x$estimate <= x$ci_upper, "estimate outside its own CI")
    chk(x$ci_level > 0 & x$ci_level < 1, "ci_level must be a proportion, e.g. 0.95")
    chk(x$comparator_type %in% c("active", "non_user"), "comparator_type must be 'active' or 'non_user'")
    chk(x$effect_measure %in% c("HR", "RR", "OR", "IRR"), "effect_measure must be HR, RR, OR or IRR")
    chk(x$brain_penetrance %in% c("high", "low", "not_applicable"), "brain_penetrance must be high, low or not_applicable")
    chk(x$drug_class %in% names(cfg$drug_classes), "drug_class not listed in config/metaanalysis.yaml drug_classes")
    for (d in paste0("robins_d", 1:7)) chk(x[[d]] %in% .robins_levels, sprintf("%s must be one of %s", d, paste(.robins_levels, collapse = "/")))
    chk(tolower(x$verified) %in% c("true", "false"), "verified must be true or false")
    chk(tolower(x$machine_rechecked) %in% c("true", "false", "discrepancy"), "machine_rechecked must be true, false or discrepancy")
    chk(tolower(x$invert_ratio) %in% c("true", "false"), "invert_ratio must be true or false")
    chk(tolower(x$include) %in% c("true", "false"), "include must be true or false")
    chk((tolower(x$include) == "true") == (x$exclusion_reason == "none"), "exclusion_reason must be 'none' exactly when include is true")
  }
  if (length(problems)) stop(paste(c(sprintf("%s failed validation:", path), problems), collapse = "\n  "))
  x
}

# Log scale and SE from a transcribed CI: SE = (log U - log L) / (2 z).
# Rows with invert_ratio = true were printed with exposure and comparator the
# other way round (e.g. Z-drug vs DORA for the DORA class): yi -> -yi; the SE
# is unchanged. The transcribed numbers themselves are never edited.
ma_to_log_scale <- function(x, tol) {
  z <- stats::qnorm(1 - (1 - x$ci_level) / 2)
  lo <- log(x$ci_lower); hi <- log(x$ci_upper)
  x$yi <- log(x$estimate)
  x$sei <- (hi - lo) / (2 * z)
  half <- (hi - lo) / 2
  off <- abs(x$yi - (lo + hi) / 2)
  # largest shift rounding to the printed decimals can cause (absent: 2 decimals)
  d <- function(col, fallback = 2L) if (is.null(x[[col]])) fallback else x[[col]]
  round_shift <- 0.5 * 10^(-d("dec_est")) / x$estimate +
    (0.5 * 10^(-d("dec_lo")) / x$ci_lower + 0.5 * 10^(-d("dec_hi")) / x$ci_upper) / 2
  bad <- which(half > 0 & off > tol * half + round_shift)
  if (length(bad)) {
    stop(sprintf("CI not symmetric on the log scale (tolerance %s of half-width) for: %s - check transcription",
                 tol, paste(if ("drug" %in% names(x)) sprintf("%s/%s", x$study_id[bad], x$drug[bad]) else x$study_id[bad],
                            collapse = ", ")))
  }
  if (any(x$sei <= 0)) stop("zero-width confidence interval in the extraction")
  # log-scale interval as transcribed (used to show a single study unchanged)
  x$lci <- lo; x$uci <- hi
  if ("invert_ratio" %in% names(x)) {
    inv <- tolower(x$invert_ratio) == "true"
    x$yi[inv] <- -x$yi[inv]
    x$lci[inv] <- -hi[inv]; x$uci[inv] <- -lo[inv]
  }
  x
}

# Random-effects pool of one group. Returns a one-row data.frame on the log scale.
ma_pool <- function(x, cfg, label) {
  k <- nrow(x)
  base <- data.frame(group = label, k = k, studies = paste(unique(x$study_id), collapse = ";"),
                     stringsAsFactors = FALSE)
  if (k < cfg$pooling$min_studies) {
    # A single study is shown with its own transcribed interval, not re-derived.
    lci <- if (is.null(x$lci)) log(x$ci_lower) else x$lci
    uci <- if (is.null(x$uci)) log(x$ci_upper) else x$uci
    return(cbind(base, est = if (k == 1) x$yi else NA_real_, ci_lower = if (k == 1) lci else NA_real_,
                 ci_upper = if (k == 1) uci else NA_real_, pi_lower = NA_real_, pi_upper = NA_real_,
                 tau2 = NA_real_, I2 = NA_real_, Q = NA_real_, Q_pval = NA_real_, est_dl = NA_real_,
                 status = sprintf("not pooled: k = %d < %d (single-study estimate shown)", k, cfg$pooling$min_studies)))
  }
  fit <- ma_rma(x$yi, x$sei, cfg)
  pr <- stats::predict(fit)
  dl <- ma_rma(x$yi, x$sei, cfg, method = "DL")
  cbind(base, est = as.numeric(fit$b), ci_lower = fit$ci.lb, ci_upper = fit$ci.ub,
        pi_lower = pr$pi.lb, pi_upper = pr$pi.ub, tau2 = fit$tau2, I2 = fit$I2 / 100,
        Q = fit$QE, Q_pval = fit$QEp, est_dl = as.numeric(dl$b), status = "pooled")
}

# Pool within each level of `by`, plus a meta-regression test for differences.
ma_subgroup <- function(x, by, cfg) {
  parts <- lapply(split(x, x[[by]]), function(g) ma_pool(g, cfg, sprintf("%s=%s", by, g[[by]][1])))
  out <- do.call(rbind, parts)
  test <- NA_real_
  if (length(unique(x[[by]])) > 1 && nrow(x) >= 3) {
    # An explicit moderator matrix (metafor adds the intercept): a formula
    # built from local variables is not evaluated in this function's scope.
    mods <- stats::model.matrix(~ factor(g), data.frame(g = x[[by]]))[, -1, drop = FALSE]
    mr <- ma_rma(x$yi, x$sei, cfg, mods = mods)
    test <- mr$QMp
  }
  out$subgroup_difference_pval <- test
  rownames(out) <- NULL
  out
}

ma_publication_bias <- function(x, cfg, label) {
  if (nrow(x) < cfg$publication_bias$min_studies) {
    return(data.frame(group = label, k = nrow(x), egger_pval = NA_real_, trimfill_k0 = NA_integer_,
                      trimfill_est = NA_real_, status = sprintf("not assessed: k = %d < %d", nrow(x), cfg$publication_bias$min_studies)))
  }
  fit <- ma_rma(x$yi, x$sei, cfg, knha = FALSE)
  eg <- metafor::regtest(fit)
  tf <- metafor::trimfill(ma_rma(x$yi, x$sei, cfg, method = "FE", knha = FALSE))
  data.frame(group = label, k = nrow(x), egger_pval = eg$pval, trimfill_k0 = tf$k0,
             trimfill_est = as.numeric(tf$b), status = "assessed")
}

ma_leave_one_out <- function(x, cfg, label) {
  if (nrow(x) < 3) return(NULL)
  fit <- ma_rma(x$yi, x$sei, cfg)
  l1 <- metafor::leave1out(fit)
  data.frame(group = label, omitted = x$study_id, est = l1$estimate, ci_lower = l1$ci.lb,
             ci_upper = l1$ci.ub, I2 = l1$I2 / 100)
}

# Positive control: the pooled anticholinergic estimate must exclude 1 on the
# known (harmful) side. Returns list(passed, detail).
ma_positive_control <- function(x, cfg) {
  pc <- cfg$positive_control
  rows <- x[x$drug_class == pc$drug_class & x$comparator_type == "non_user", , drop = FALSE]
  if (nrow(rows) == 0) return(list(passed = FALSE, detail = sprintf("no '%s' rows in the extraction", pc$drug_class)))
  p <- ma_pool(rows, cfg, "positive_control")
  if (p$status != "pooled") return(list(passed = FALSE, detail = p$status))
  ok <- if (pc$expected_direction == "harmful") p$ci_lower > 0 else p$ci_upper < 0
  list(passed = ok, pool = p,
       detail = sprintf("%s: pooled ratio %.2f [%.2f, %.2f], k = %d", pc$drug_class, exp(p$est),
                        exp(p$ci_lower), exp(p$ci_upper), p$k))
}

# Sign agreement between two log-scale estimates/intervals.
.agree <- function(a, a_lo, a_hi, expected_sign) {
  if (is.na(a) || is.na(expected_sign)) return(if (is.na(a)) "missing" else "no prediction")
  if (a_lo > 0 || a_hi < 0) { if (sign(a) == expected_sign) "agree" else "DISAGREE" } else "uncertain (CI includes null)"
}

# Triangulation: per MR target, the pooled observational estimate, the MR
# estimate and the mechanistic prediction, all expressed as the effect of
# LOWER target activity on log risk. Missing pieces stay NA (shown as gaps).
ma_triangulate <- function(pooled_by_class, mr_summary, mr_targets, cfg) {
  rows <- list()
  for (gene in names(mr_targets)) {
    pred <- mr_targets[[gene]]$predicted_lower_activity
    exp_sign <- switch(pred, protective = -1, harmful = 1, exploratory = NA)
    classes <- names(Filter(function(d) identical(d$target, gene), cfg$drug_classes))
    obs <- pooled_by_class[pooled_by_class$drug_class %in% classes & !is.na(pooled_by_class$est), , drop = FALSE]
    if (nrow(obs) > 1) stop(sprintf("%s: more than one pooled drug class maps to the target; pool or choose explicitly", gene))
    if (nrow(obs) == 1) {
      flip <- if (cfg$drug_classes[[obs$drug_class]]$action == "activates") -1 else 1
      o_est <- obs$est * flip
      o <- sort(c(obs$ci_lower, obs$ci_upper) * flip)
    } else {
      o_est <- NA_real_; o <- c(NA_real_, NA_real_)
    }
    m <- mr_summary[mr_summary$gene == gene & mr_summary$analysis %in% "main" & mr_summary$is_primary %in% TRUE, , drop = FALSE]
    if (nrow(m) > 1) m <- m[m$exposure_source == m$exposure_source[1], , drop = FALSE]  # first source listed
    mr <- if (nrow(m) == 1) c(-m$b, sort(-c(m$ci_lower, m$ci_upper))) else c(NA, NA, NA)
    rows[[gene]] <- data.frame(
      target = gene, prediction_lower_activity = pred,
      obs_est = o_est, obs_ci_lower = o[1], obs_ci_upper = o[2],
      obs_k = if (nrow(obs) == 1) obs$k else 0L, obs_drug_class = if (nrow(obs) == 1) obs$drug_class else NA,
      mr_est = mr[1], mr_ci_lower = mr[2], mr_ci_upper = mr[3],
      mr_source = if (nrow(m) == 1) m$exposure_source else NA,
      obs_vs_prediction = .agree(o_est, o[1], o[2], exp_sign),
      mr_vs_prediction = .agree(mr[1], mr[2], mr[3], exp_sign),
      obs_vs_mr = if (is.na(o_est) || is.na(mr[1])) "missing" else if (sign(o_est) == sign(mr[1])) "same sign" else "OPPOSITE SIGN",
      stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  out
}
