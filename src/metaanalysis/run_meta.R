# Module 5 entry point. Usage: Rscript src/metaanalysis/run_meta.R  (or `make metaanalysis`)
#
# Refuses to start unless the Module 5 validation passed on the current code
# and prereg/metaanalysis.md is frozen. The anticholinergic positive control
# is checked FIRST: if the pipeline cannot recover that known association it
# is broken, and the run stops before any noradrenergic result is produced.

options(warn = 2)
source("src/common/provenance.R")
source("src/common/io.R")
source("src/metaanalysis/meta.R")

root <- lc_repo_root()
cfg_path <- file.path(root, "config", "metaanalysis.yaml")
cfg <- lc_read_yaml(cfg_path)
lc_require_validated("metaanalysis", root)
prereg <- lc_require_preregistration("metaanalysis", root)

lc_with_run("metaanalysis", cfg_path, seed = cfg$seed, packages = c("metafor", "yaml", "digest"),
            preregistration = prereg, body = function(run) {
  tables <- file.path(root, "results", "tables")
  dir.create(tables, recursive = TRUE, showWarnings = FALSE)
  write_out <- function(df, name) { p <- file.path(tables, name); utils::write.csv(df, p, row.names = FALSE); lc_record_output(run, p) }

  x <- ma_load_extraction(lc_analysis_input(file.path(root, cfg$extraction$file), "researcher, hand extraction", root), cfg)
  x <- ma_to_log_scale(x, cfg$extraction$log_symmetry_tolerance)
  lc_log(run, sprintf("extraction: %d rows from %d studies", nrow(x), length(unique(x$study_id))))

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

  by_class <- do.call(rbind, lapply(split(ratio, list(ratio$drug_class, ratio$comparator_type), drop = TRUE), function(g) {
    p <- ma_pool(g, cfg, sprintf("%s vs %s", g$drug_class[1], g$comparator_type[1]))
    cbind(drug_class = g$drug_class[1], comparator_type = g$comparator_type[1], p)
  }))
  write_out(by_class, "meta_pooled_by_class.csv")
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

  # Triangulation with Module 1 (MR); missing MR results stay as gaps.
  mr_path <- file.path(tables, "mr_summary.csv")
  mr <- if (file.exists(mr_path)) utils::read.csv(mr_path, stringsAsFactors = FALSE) else {
    lc_log(run, "results/tables/mr_summary.csv absent: MR column of the triangulation is empty", "WARNING")
    data.frame(gene = character(), analysis = character(), is_primary = logical(), b = numeric(),
               ci_lower = numeric(), ci_upper = numeric(), exposure_source = character())
  }
  tri_in <- by_class[by_class$comparator_type == cfg$triangulation$comparator_type, , drop = FALSE]
  mr_targets <- lc_read_yaml(file.path(root, "config", "mr_targets.yaml"))$targets
  write_out(ma_triangulate(tri_in, mr, mr_targets, cfg), "triangulation.csv")
})
