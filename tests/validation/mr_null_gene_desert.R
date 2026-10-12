# Null test: random gene-desert variant sets as pseudo-instruments for a
# synthetic exposure, against the real AD outcome. See the target YAML for the
# design. Returns list(passed, detail, table).

validate_mr_null_gene_desert <- function(run, settings, root = lc_repo_root()) {
  target <- lc_load_target("mr_null_gene_desert", root)
  s <- target$settings
  mr_cfg <- lc_read_yaml(file.path(root, "config", "mr_targets.yaml"))
  if (!identical(target$inputs$genome_build, mr_cfg$genome_build)) {
    stop(sprintf("null test build %s != Module 1 build %s", target$inputs$genome_build, mr_cfg$genome_build))
  }
  deserts <- utils::read.csv(lc_analysis_input(file.path(root, target$inputs$desert_variants),
                                               target$inputs$provided_by, root), stringsAsFactors = FALSE)
  if (!"SNP" %in% names(deserts)) stop("desert variant file needs a SNP column")
  out <- mr_read_sumstats(mr_cfg$outcome$file, mr_cfg$outcome$columns, mr_cfg$genome_build,
                          "researcher download (AD outcome)", NULL, NULL, root)
  avail <- out$data[out$data$SNP %in% deserts$SNP, , drop = FALSE]
  lc_log(run, sprintf("null test: %d of %d desert variants present in the outcome", nrow(avail), nrow(deserts)))
  if (nrow(avail) < s$variants_per_set * 10) stop("too few desert variants in the outcome for a meaningful null test")

  z <- stats::qnorm(1 - (1 - settings$ci_level) / 2)
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(s$seed)
  p <- vapply(seq_len(s$n_sets), function(i) {
    o <- avail[sample.int(nrow(avail), s$variants_per_set), ]
    h <- data.frame(beta_exp = stats::runif(nrow(o), 0.05, 0.2), se_exp = 0.01,
                    beta_out = o$beta, se_out = o$se)
    mr_ivw(h, z)$pval
  }, numeric(1))
  rate <- mean(p < s$alpha)
  limit <- stats::qbinom(s$pass_quantile, s$n_sets, s$alpha) / s$n_sets
  list(passed = rate <= limit,
       detail = sprintf("IVW p < %g in %.1f%% of %d null sets (limit %.1f%%)", s$alpha, 100 * rate, s$n_sets, 100 * limit),
       table = data.frame(target = target$id, n_sets = s$n_sets, variants_per_set = s$variants_per_set,
                          alpha = s$alpha, false_positive_rate = rate, pass_limit = limit,
                          ks_uniform_pval = stats::ks.test(p, "punif")$p.value))
}
