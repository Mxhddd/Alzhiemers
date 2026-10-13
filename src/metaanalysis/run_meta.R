# Module 5 entry point.
#   Rscript src/metaanalysis/run_meta.R                      (real extraction; gated)
#   Rscript src/metaanalysis/run_meta.R --dry-run-simulated  (simulated; no gates)
#
# Real mode refuses to start unless the Module 5 validation passed on the
# current code and prereg/metaanalysis.md is frozen. The extraction is
# data/extracted/meta_extraction.csv (Claude-extracted from published papers,
# verified=false until a human checks it: every table is then stamped
# UNVERIFIED). The dry run builds a simulated extraction in memory
# (sim/meta_sim.R), writes to results/dryrun/ and stamps everything SIMULATED.
# The anticholinergic positive control is checked FIRST: if the pipeline
# cannot recover that known association it is broken, and the run stops
# before any noradrenergic result is produced.

options(warn = 2)
source("src/common/provenance.R")
source("src/common/io.R")
source("src/metaanalysis/meta.R")

dry <- "--dry-run-simulated" %in% commandArgs(trailingOnly = TRUE)
root <- lc_repo_root()
cfg_path <- file.path(root, "config", "metaanalysis.yaml")
cfg <- lc_read_yaml(cfg_path)
prereg <- NULL
if (!dry) {
  lc_require_validated("metaanalysis", root)
  prereg <- lc_require_preregistration("metaanalysis", root)
}

lc_with_run(if (dry) "metaanalysis_dryrun" else "metaanalysis", cfg_path, seed = cfg$seed,
            packages = c("metafor", "yaml", "digest"), preregistration = prereg, body = function(run) {
  tables <- file.path(root, "results", if (dry) "dryrun/tables" else "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  write_out <- function(df, name) lc_write_table(run, df, file.path(tables, name))

  if (dry) {
    lc_mark_simulated(run)
    source(file.path(root, "sim", "meta_sim.R"))
    p <- tempfile(fileext = ".csv")
    utils::write.csv(meta_simulate_extraction(cfg, cfg$seed), p, row.names = FALSE)
    x <- ma_load_extraction(p, cfg)
  } else {
    path <- lc_analysis_input(file.path(root, cfg$extraction$file), "Claude extraction from published papers", root)
    lc_load_extracted(path, "Claude extraction from published papers", run, root)   # provenance + UNVERIFIED stamp
    x <- ma_load_extraction(path, cfg)
  }
  # Rows the extraction itself marks include=false (with a reason) are counted
  # and logged as drops, never silently skipped.
  x_all <- x                                   # kept for the E1 sensitivity analysis
  excl <- tolower(x$include) != "true"
  for (i in which(excl)) lc_drop(run, "extraction", sprintf("%s (%s)", x$study_id[i], x$drug_class[i]), x$exclusion_reason[i])
  x <- x[!excl, , drop = FALSE]
  x <- ma_to_log_scale(x, cfg$extraction$log_symmetry_tolerance)
  lc_log(run, sprintf("extraction: %d included rows from %d studies (%d excluded rows logged)", nrow(x),
                      length(unique(x$study_id)), sum(excl)))

  # Effect measures that cannot be combined are analysed separately, never dropped.
  x$measure_group <- ifelse(x$effect_measure %in% cfg$pooling$combinable_measures, "ratio", x$effect_measure)
  for (g in setdiff(unique(x$measure_group), "ratio")) {
    lc_log(run, sprintf("%d rows with effect measure %s pooled separately from %s",
                        sum(x$measure_group == g), g, paste(cfg$pooling$combinable_measures, collapse = "/")), "WARNING")
  }

  ratio <- x[x$measure_group == "ratio", , drop = FALSE]
  pc <- ma_positive_control(ratio, cfg)
  lc_log(run, sprintf("positive control (%s): %s - %s", cfg$positive_control$drug_class,
                      if (pc$passed) "RECOVERED" else "NOT RECOVERED", pc$detail), if (pc$passed) "INFO" else "ERROR")
  if (!pc$passed) stop("anticholinergic positive control not recovered: the pipeline is broken or the extraction is wrong; no results produced")

  # Active-comparator rows are pooled within their contrast only (pooling
  # "vs hydrophilic beta-blocker" with "vs supplements" would mix questions).
  ratio$pool_contrast <- ifelse(ratio$comparator_type == "active", ratio$contrast, "none")
  by_class <- do.call(rbind, lapply(split(ratio, list(ratio$drug_class, ratio$comparator_type, ratio$pool_contrast), drop = TRUE),
                                    function(g) {
    p <- ma_pool(g, cfg, sprintf("%s vs %s%s", g$drug_class[1], g$comparator_type[1],
                                 if (g$comparator_type[1] == "active") sprintf(" [%s]", g$pool_contrast[1]) else ""))
    cbind(drug_class = g$drug_class[1], comparator_type = g$comparator_type[1], contrast = g$pool_contrast[1], p)
  }))
  write_out(by_class, "meta_pooled_by_class.csv")

  # Pre-registered contrast verdicts (prereg §2) and their sensitivity analyses (§8).
  verdicts <- ma_contrast_verdicts(ratio, cfg, "primary")
  sens <- cfg$sensitivity
  if (isTRUE(sens$dl_estimator)) {
    c2 <- cfg; c2$pooling$method <- "DL"; verdicts <- rbind(verdicts, ma_contrast_verdicts(ratio, c2, "DerSimonian-Laird"))
  }
  if (isTRUE(sens$no_knapp_hartung)) {
    c2 <- cfg; c2$pooling$knapp_hartung <- FALSE; verdicts <- rbind(verdicts, ma_contrast_verdicts(ratio, c2, "no Knapp-Hartung"))
  }
  if (isTRUE(sens$combine_or)) {
    verdicts <- rbind(verdicts, ma_contrast_verdicts(x, cfg, "ORs pooled with HR/RR/IRR"))
  }
  if (isTRUE(sens$include_e1_rows)) {
    e1 <- x_all[tolower(x_all$include) != "true" & startsWith(x_all$exclusion_reason, "E1"), , drop = FALSE]
    if (nrow(e1)) {
      e1 <- ma_to_log_scale(e1, cfg$extraction$log_symmetry_tolerance)
      e1$measure_group <- ifelse(e1$effect_measure %in% cfg$pooling$combinable_measures, "ratio", e1$effect_measure)
      common <- intersect(names(ratio), names(e1))
      verdicts <- rbind(verdicts, ma_contrast_verdicts(rbind(ratio[, common], e1[e1$measure_group == "ratio", common]), cfg,
                                                       "E1 rows (Parkinson's populations) included"))
    }
  }
  write_out(verdicts, "meta_contrast_verdicts.csv")
  for (i in which(verdicts$analysis == "primary")) {
    lc_log(run, sprintf("PRE-REGISTERED CONTRAST %s: %s (k = %d)", verdicts$contrast[i], verdicts$verdict[i], verdicts$k[i]))
  }
  write_out(ma_subgroup(ratio[ratio$contrast %in% cfg$contrasts, , drop = FALSE], "contrast", cfg), "meta_contrasts.csv")
  write_out(ma_subgroup(ratio, "comparator_type", cfg), "meta_by_comparator_type.csv")
  pen <- ratio[ratio$brain_penetrance != "not_applicable", , drop = FALSE]
  if (nrow(pen)) write_out(ma_subgroup(pen, "brain_penetrance", cfg), "meta_by_brain_penetrance.csv")
  others <- x[x$measure_group != "ratio", , drop = FALSE]
  if (nrow(others)) write_out(do.call(rbind, lapply(split(others, list(others$drug_class, others$effect_measure), drop = TRUE),
                                                    function(g) ma_pool(g, cfg, sprintf("%s (%s)", g$drug_class[1], g$effect_measure[1])))),
                              "meta_pooled_other_measures.csv")

  groups <- split(ratio, ratio$drug_class)
  write_out(do.call(rbind, lapply(names(groups), function(n) ma_publication_bias(groups[[n]], cfg, n))), "meta_publication_bias.csv")
  loo <- do.call(rbind, lapply(names(groups), function(n) ma_leave_one_out(groups[[n]], cfg, n)))
  if (!is.null(loo)) write_out(loo, "meta_leave_one_out.csv")
  write_out(x, "meta_extraction_log_scale.csv")
  lc_write_drop_report(run, file.path(root, "results", if (dry) "dryrun/logs" else "logs", "meta_dropped.md"))

  # Triangulation with Module 1 (MR); missing MR results stay as gaps.
  mr_path <- file.path(root, "results", "tables", "mr_summary.csv")
  mr <- if (file.exists(mr_path)) utils::read.csv(mr_path, stringsAsFactors = FALSE) else {
    lc_log(run, "results/tables/mr_summary.csv absent: MR column of the triangulation is empty", "WARNING")
    data.frame(gene = character(), analysis = character(), is_primary = logical(), b = numeric(),
               ci_lower = numeric(), ci_upper = numeric(), exposure_source = character())
  }
  tri_in <- by_class[by_class$comparator_type == cfg$triangulation$comparator_type &
                       (by_class$comparator_type != "active" | by_class$contrast %in% unlist(cfg$contrasts)), , drop = FALSE]
  mr_targets <- lc_read_yaml(file.path(root, "config", "mr_targets.yaml"))$targets
  write_out(ma_triangulate(tri_in, mr, mr_targets, cfg), "triangulation.csv")
})
