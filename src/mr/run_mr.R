# Module 1 entry point: drug-target MR for every pre-registered target.
# Usage: Rscript src/mr/run_mr.R   (or `make mr`)
#
# Refuses to start unless (1) the MR validation passed on the current code and
# (2) prereg/mr.md is frozen (python -m src.common.prereg check mr).
# Outputs: results/tables/mr_summary.csv, results/logs/mr_dropped.md, and the
# forest plot (drawn by src/mr/forest_plot.py with the shared theme).

options(warn = 2)  # fail loudly: every warning is an error
source("src/common/provenance.R")
source("src/common/io.R")
for (f in c("harmonise", "estimators", "crosscheck", "regions", "sumstats", "instruments", "coloc", "steiger", "targets")) {
  source(file.path("src", "mr", paste0(f, ".R")))
}

root <- lc_repo_root()
cfg_path <- file.path(root, "config", "mr_targets.yaml")
cfg <- lc_read_yaml(cfg_path)
mr_settings <- lc_read_yaml(file.path(root, "config", "mr_methods.yaml"))

# ---- Gates -------------------------------------------------------------------
lc_require_validated("mr", root)
prereg <- lc_require_preregistration("mr", root)

# ---- Researcher-supplied fields that the main analysis cannot run without ----
need <- list("outcome.file" = cfg$outcome$file, "outcome.columns" = cfg$outcome$columns,
             "outcome.n_cases" = cfg$outcome$n_cases, "outcome.n_controls" = cfg$outcome$n_controls,
             "instruments.ld_reference_bfile" = cfg$instruments$ld_reference_bfile)
missing <- names(need)[vapply(need, is.null, logical(1))]
if (length(missing)) stop(sprintf("config/mr_targets.yaml: the researcher must supply %s", paste(missing, collapse = ", ")))
mr_check_targets(cfg$targets)
for (s in names(cfg$exposures)) {
  b <- cfg$exposures[[s]]$build
  if (!is.null(cfg$exposures[[s]]$file) && !identical(b, cfg$genome_build)) {
    stop(sprintf("exposure %s is on build %s but the outcome is %s: lift over first", s, b, cfg$genome_build))
  }
}

io <- list(
  r2_fun = mr_plink_r2_fun(file.path(root, cfg$instruments$ld_reference_bfile)),
  exposure = function(source_name, region, gene) {
    e <- cfg$exposures[[source_name]]
    mr_read_sumstats(e$file, e$columns, e$build, sprintf("researcher download (%s)", source_name), region, gene, root)
  },
  positive_control = function(gene, pc, region) {
    if (is.null(pc$file)) stop(sprintf("positive control %s for %s has no file (researcher to supply)", pc$trait, gene))
    mr_read_sumstats(pc$file, pc$columns, cfg$genome_build, sprintf("researcher download (%s GWAS)", pc$trait), region, NULL, root)
  },
  outcome = function(which, region) {
    spec <- switch(which,
      main = list(file = cfg$outcome$file, columns = cfg$outcome$columns),
      proxy_free = list(file = cfg$outcome$proxy_free_file, columns = cfg$outcome$columns),
      replication = list(file = cfg$replication$file, columns = cfg$replication$columns))
    if (is.null(spec$file)) { if (which == "main") stop("no outcome file"); return(NULL) }
    mr_read_sumstats(spec$file, spec$columns, cfg$genome_build, sprintf("researcher download (%s)", which), region, NULL, root)
  }
)

lc_with_run("mr", cfg_path, seed = cfg$seed,
            packages = c("data.table", "biomaRt", "yaml", "digest", "TwoSampleMR"),
            preregistration = prereg, body = function(run) {
  drop <- function(stage, item, reason) lc_drop(run, stage = stage, item = item, reason = reason)
  logf <- function(msg) lc_log(run, msg)
  parts <- list()
  cc_parts <- list()
  # Primary window (per target) and the pre-registered sensitivity window.
  for (gene in names(cfg$targets)) {
    tcfg <- cfg$targets[[gene]]
    for (wkb in c(tcfg$window_kb, cfg$window_kb_sensitivity)) {
      region <- mr_gene_region(gene, wkb, run = run, root = root)
      for (s in names(cfg$exposures)) {
        if (is.null(cfg$exposures[[s]]$file)) {
          drop("exposure_source", sprintf("%s [%s, +/-%s kb]", gene, s, wkb), "exposure file not supplied")
          parts[[length(parts) + 1]] <- data.frame(gene = gene, exposure_source = s, window_kb = wkb,
            predicted_lower_activity = tcfg$predicted_lower_activity, method = NA_character_,
            status = "excluded: exposure file not supplied", stringsAsFactors = FALSE)
          next
        }
        res <- mr_analyse_target(gene, tcfg, s, region, io, cfg, mr_settings, drop, logf)
        res$rows$window_kb <- wkb
        parts[[length(parts) + 1]] <- res$rows
        if (!is.null(res$crosscheck)) {
          res$crosscheck$window_kb <- wkb
          cc_parts[[length(cc_parts) + 1]] <- res$crosscheck
        }
      }
    }
  }
  summary <- mr_bind_rows(parts)
  # Primary source per target and window (prereg section 5.2): the first source
  # in the hierarchy whose instruments passed the positive control (rows that
  # were not excluded). Blood-only rows are flagged as option (a).
  passed_gate <- !grepl("^excluded", summary$status)
  gate <- unique(data.frame(gene = summary$gene, window_kb = summary$window_kb,
                            exposure_source = summary$exposure_source, passed = passed_gate))
  sel <- mr_select_source(gate, cfg$source_selection$hierarchy)
  key <- paste(summary$gene, summary$window_kb)
  summary$primary_source <- sel$primary_source[match(key, paste(sel$gene, sel$window_kb))]
  summary$is_primary_source <- !is.na(summary$primary_source) & summary$exposure_source == summary$primary_source
  summary$blood_only_option_a <- summary$exposure_source == cfg$source_selection$blood_only_source
  out <- file.path(root, "results", "tables", "mr_summary.csv")
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(summary, out, row.names = FALSE)
  lc_record_output(run, out)
  # Pre-registered cross-check against TwoSampleMR: every comparison, agreeing or not
  cc_out <- file.path(root, "results", "tables", "mr_crosscheck_twosamplemr.csv")
  cc_all <- if (length(cc_parts)) do.call(rbind, cc_parts) else data.frame()
  utils::write.csv(cc_all, cc_out, row.names = FALSE)
  lc_record_output(run, cc_out)
  if (nrow(cc_all) && any(!cc_all$agree)) {
    lc_log(run, level = "WARNING", sprintf("%d base-R vs TwoSampleMR comparison(s) disagree; the TwoSampleMR value was used (%s)",
                         sum(!cc_all$agree), cc_out))
  }
  lc_write_drop_report(run, file.path(root, "results", "logs", "mr_dropped.md"))
})
