# Validation harness, R side. A validation target is a YAML file in
# tests/validation/targets/ holding a published result transcribed by the
# researcher. The harness refuses incomplete targets, reads real inputs only
# through lc_analysis_input(), and compares our estimate with the published
# confidence interval on the log scale.

# Required fields by target kind.
lc_target_required <- list(
  published_estimate = list(
    published = c("citation", "doi", "locator", "method", "scale", "reported_as",
                  "estimate", "ci_lower", "ci_upper", "ci_level", "transcribed_by", "transcribed_on"),
    conversion = c("factor", "reason"),
    inputs = c("genome_build", "exposure", "outcome", "provided_by")),
  null_test = list(
    settings = c("seed", "n_sets", "variants_per_set", "alpha", "pass_quantile"),
    inputs = c("genome_build", "desert_variants", "provided_by")),
  # The pooled estimate, its CI and the component studies are in the
  # extracted data file named by inputs.studies (role ma_pooled / ma_component,
  # each row with source, locator and verbatim snippet), not typed in the YAML.
  meta_published = list(
    published = c("citation", "doi", "locator", "method", "knapp_hartung", "reported_as",
                  "ci_level", "transcribed_by", "transcribed_on"),
    inputs = c("studies", "provided_by")),
  positive_control = list(),
  hand_arithmetic = list(
    inputs = c("yi", "vi", "provided_by"),
    expected = c("fe_estimate", "fe_se", "q", "tau2_dl", "re_dl_estimate", "re_dl_se", "tolerance")),
  simulation = list(
    settings = c("seed", "n_reps", "k", "true_log_effect", "se_min", "se_max", "max_median_i2"))
)

lc_load_target <- function(id, root = lc_repo_root()) {
  path <- file.path(root, "tests", "validation", "targets", paste0(id, ".yaml"))
  if (!file.exists(path)) stop(sprintf("validation target %s not found", path))
  t <- lc_read_yaml(path)
  kind <- t$kind
  if (is.null(kind) || !kind %in% names(lc_target_required)) {
    stop(sprintf("target '%s': kind must be one of %s", id, paste(names(lc_target_required), collapse = ", ")))
  }
  req <- lc_target_required[[kind]]
  missing <- character()
  for (section in names(req)) {
    for (field in req[[section]]) {
      v <- t[[section]][[field]]
      if (is.null(v) || (length(v) == 1 && is.na(v)) || identical(v, "")) {
        missing <- c(missing, paste(section, field, sep = "."))
      }
    }
  }
  if (length(missing)) {
    stop(sprintf(paste0("validation target '%s' is incomplete - missing: %s (in %s). Published values ",
                        "must come from the source paper with a locator (data/extracted/), never be guessed."),
                 id, paste(missing, collapse = ", "), lc_relpath(path, root)))
  }
  if (kind != "published_estimate") return(t)
  p <- t$published
  if (!p$reported_as %in% c("ratio", "log")) stop("published.reported_as must be 'ratio' or 'log'")
  if (!(p$ci_lower <= p$estimate && p$estimate <= p$ci_upper)) {
    stop(sprintf("target '%s': published estimate is not inside its own CI - check the transcription", id))
  }
  if (p$reported_as == "ratio" && any(c(p$estimate, p$ci_lower, p$ci_upper) <= 0)) {
    stop(sprintf("target '%s': ratio-scale values must be positive", id))
  }
  t
}

lc_read_sumstats <- function(path, provided_by, genome_build, root = lc_repo_root()) {
  abs <- lc_analysis_input(path, provided_by, root = root)
  x <- utils::read.csv(abs, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  attr(x, "genome_build") <- genome_build
  mr_check_sumstats(x, path)
  x
}

# Compare one estimate row (method, b, se) with the published interval.
lc_compare_to_target <- function(row, target, our_ci_level) {
  p <- target$published
  to_log <- function(v) if (p$reported_as == "ratio") log(v) else v
  f <- target$conversion$factor
  ours <- row$b * f
  ours_ci <- sort(c(row$ci_lower, row$ci_upper) * f)
  lo <- to_log(p$ci_lower); hi <- to_log(p$ci_upper)
  data.frame(
    target = target$id, method = row$method, nsnp = row$nsnp,
    our_estimate_log = ours, our_ci_lower_log = ours_ci[1], our_ci_upper_log = ours_ci[2],
    our_ci_level = our_ci_level,
    published_estimate_log = to_log(p$estimate), published_ci_lower_log = lo,
    published_ci_upper_log = hi, published_ci_level = p$ci_level,
    within_published_ci = lo <= ours && ours <= hi,
    published_within_our_ci = ours_ci[1] <= to_log(p$estimate) && to_log(p$estimate) <= ours_ci[2],
    stringsAsFactors = FALSE
  )
}
