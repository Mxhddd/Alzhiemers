# Module 1, positive-control step ONLY: for every pre-registered target, build
# the cis instruments from each configured exposure and run the MR on the
# target's positive-control trait (config/mr_targets.yaml). This is the gate of
# prereg/mr.md §3, evaluated on its own.
#
# Usage: Rscript src/mr/run_positive_controls.R
#
# It never reads an Alzheimer's outcome: there is no outcome reader here, and
# the script stops if any file it is about to read lies under
# data/raw/mr/quarantine_outcome/. It does not need the frozen pre-registration
# (no AD data are involved) and it does not change what run_mr.R does: the gate
# is re-evaluated there, with the same functions, before the outcome is read.
#
# Outputs:
#   results/tables/mr_positive_controls.csv            one row per target x exposure x window
#   results/tables/mr_positive_controls_crosscheck.csv base-R vs TwoSampleMR, every comparison
#   results/logs/mr_positive_controls_dropped.md       every dropped variant with its reason

options(warn = 2)  # fail loudly: every warning is an error
source("src/common/provenance.R")
source("src/common/io.R")
for (f in c("harmonise", "estimators", "crosscheck", "regions", "sumstats", "instruments", "targets")) {
  source(file.path("src", "mr", paste0(f, ".R")))
}

root <- lc_repo_root()
cfg_path <- file.path(root, "config", "mr_targets.yaml")
cfg <- lc_read_yaml(cfg_path)
mr_settings <- lc_read_yaml(file.path(root, "config", "mr_methods.yaml"))
mr_check_targets(cfg$targets)
if (is.null(cfg$instruments$ld_reference_bfile)) stop("config/mr_targets.yaml: instruments.ld_reference_bfile is null")

quarantine <- normalizePath(file.path(root, "data", "raw", "mr", "quarantine_outcome"), mustWork = FALSE)
assert_not_outcome <- function(path) {
  abs <- normalizePath(file.path(root, path), mustWork = FALSE)
  if (startsWith(abs, paste0(quarantine, "/"))) {
    stop(sprintf("%s is a quarantined Alzheimer's outcome file; the positive-control step must not read it", path))
  }
  path
}

# Primary window (per target) and the pre-registered sensitivity window.
windows <- function(tcfg) c(primary = tcfg$window_kb, sensitivity = cfg$window_kb_sensitivity)

lc_with_run("mr_positive_controls", cfg_path, seed = cfg$seed,
            packages = c("data.table", "biomaRt", "yaml", "digest", "TwoSampleMR"), body = function(run) {
  drop <- function(stage, item, reason) lc_drop(run, stage = stage, item = item, reason = reason)
  r2_fun <- mr_plink_r2_fun(file.path(root, cfg$instruments$ld_reference_bfile))
  rows <- list()
  cc <- list()
  add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

  for (gene in names(cfg$targets)) {
    tcfg <- cfg$targets[[gene]]
    pc <- tcfg$positive_control
    for (analysis in names(windows(tcfg))) {
      wkb <- windows(tcfg)[[analysis]]
      base <- list(gene = gene, analysis = analysis, window_kb = wkb,
                   trait = if (is.null(pc)) NA_character_ else pc$trait,
                   expected_sign = if (is.null(pc)) NA_integer_ else as.integer(pc$expected_sign))
      if (is.null(pc)) {
        drop("positive_control", sprintf("%s [%s]", gene, analysis), "no positive control defined: untestable")
        add(base, exposure_source = NA_character_, status = "untestable: no positive control defined (pre-registered)")
        next
      }
      region <- mr_gene_region(gene, wkb, run = run, root = root)
      for (s in names(cfg$exposures)) {
        e <- cfg$exposures[[s]]
        b <- c(base, exposure_source = s)
        if (is.null(e$file)) {
          add(b, status = "not run: exposure file not supplied")
          next
        }
        exp_raw <- mr_read_sumstats(assert_not_outcome(e$file), e$columns, e$build,
                                    sprintf("Module 1 data step (%s)", s), region, gene, root)
        for (i in seq_len(nrow(exp_raw$dropped))) drop(sprintf("exposure_read_%s", analysis), sprintf("%s %s", gene, exp_raw$dropped$SNP[i]), exp_raw$dropped$reason[i])
        if (nrow(exp_raw$data) == 0) {
          add(b, status = sprintf("excluded: no %s variants for %s in the cis window (gene not in the exposure file)", s, gene))
          next
        }
        inst <- mr_build_instruments(exp_raw$data, cfg$instruments, r2_fun)
        for (i in seq_len(nrow(inst$dropped))) drop(sprintf("instruments_%s", analysis), sprintf("%s %s", gene, inst$dropped$SNP[i]), inst$dropped$reason[i])
        if (nrow(inst$data) == 0) {
          add(b, n_cis_variants = nrow(exp_raw$data),
              status = "excluded: no instrument passed p, MAF, F and clumping filters")
          next
        }
        lc_log(run, sprintf("%s [%s, %s, +/-%s kb]: %d instrument(s): %s", gene, s, analysis, wkb,
                            nrow(inst$data), paste(inst$data$SNP, collapse = ", ")))
        if (is.null(pc$file)) {
          add(b, n_cis_variants = nrow(exp_raw$data), n_instruments = nrow(inst$data),
              status = "not run: positive-control summary statistics not supplied")
          next
        }
        pc_ss <- mr_read_sumstats(assert_not_outcome(pc$file), pc$columns, cfg$genome_build,
                                  sprintf("Module 1 data step (%s GWAS)", pc$trait), region, NULL, root)
        for (i in seq_len(nrow(pc_ss$dropped))) drop(sprintf("control_read_%s", analysis), sprintf("%s %s", gene, pc_ss$dropped$SNP[i]), pc_ss$dropped$reason[i])
        res <- mr_positive_control(inst$data, pc_ss$data, pc, mr_settings, cfg$positive_control$p_max)
        hz <- res$harmonise
        for (i in seq_len(nrow(hz$dropped))) {
          if (hz$dropped$SNP[i] %in% inst$data$SNP) drop(sprintf("harmonise_%s", analysis), sprintf("%s %s", gene, hz$dropped$SNP[i]), hz$dropped$reason[i])
        }
        est <- res$estimate
        lc_log(run, sprintf("%s [%s, %s] positive control %s: %s", gene, s, analysis, pc$trait, res$reason))
        add(b, n_cis_variants = nrow(exp_raw$data), n_instruments = nrow(inst$data),
            instruments = paste(inst$data$SNP, collapse = ";"),
            f_min = min(inst$data$f_stat), f_median = stats::median(inst$data$f_stat),
            n_harmonised = nrow(hz$data),
            harmonise_actions = paste(names(hz$counts), hz$counts, sep = "=", collapse = ","),
            method = if (is.null(est)) NA_character_ else est$method,
            b = if (is.null(est)) NA_real_ else est$b,
            se = if (is.null(est)) NA_real_ else est$se,
            ci_lower = if (is.null(est)) NA_real_ else est$ci_lower,
            ci_upper = if (is.null(est)) NA_real_ else est$ci_upper,
            pval = if (is.null(est)) NA_real_ else est$pval,
            passed = res$passed,
            status = if (res$passed) "passed" else res$reason)
        if (!is.null(res$crosscheck)) cc[[length(cc) + 1]] <- data.frame(gene = gene, analysis = analysis, exposure_source = s, res$crosscheck)
      }
    }
  }
  out <- mr_bind_rows(rows)
  path <- file.path(root, "results", "tables", "mr_positive_controls.csv")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, path, row.names = FALSE)
  lc_record_output(run, path)
  cc_all <- if (length(cc)) do.call(rbind, cc) else data.frame()
  cc_path <- file.path(root, "results", "tables", "mr_positive_controls_crosscheck.csv")
  utils::write.csv(cc_all, cc_path, row.names = FALSE)
  lc_record_output(run, cc_path)
  if (nrow(cc_all) && any(!cc_all$agree)) {
    lc_log(run, level = "WARNING", sprintf("%d base-R vs TwoSampleMR comparison(s) disagree (%s)", sum(!cc_all$agree), cc_path))
  }
  lc_write_drop_report(run, file.path(root, "results", "logs", "mr_positive_controls_dropped.md"))
})
