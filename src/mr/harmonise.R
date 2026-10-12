# Allele harmonisation for two-sample MR.
#
# Puts exposure and outcome effects on the same effect allele. Every variant
# ends in exactly one of these actions (kept) or reasons (dropped):
#
#   kept:    direct               same alleles, same orientation
#            swapped              alleles reversed -> outcome beta negated, EAF 1-EAF
#            strand_flip          outcome on the other strand (A<->T, C<->G)
#            strand_flip_swapped  both of the above
#            palindrome_direct    A/T or C/G variant, EAFs agree on orientation
#            palindrome_flipped   A/T or C/G variant, EAFs say outcome is flipped
#   dropped: absent from outcome / absent from exposure
#            alleles incompatible
#            palindromic, EAF missing
#            palindromic, MAF above threshold (strand ambiguous)
#            non-SNV alleles that do not match directly or swapped
#
# Inputs must carry attr(x, "genome_build"); a mismatch stops the run because
# a build mismatch produces wrong instruments with no other warning.

.mr_required_cols <- c("SNP", "effect_allele", "other_allele", "beta", "se", "eaf")

mr_check_sumstats <- function(x, what) {
  if (!is.data.frame(x)) stop(sprintf("%s must be a data.frame", what))
  missing <- setdiff(.mr_required_cols, names(x))
  if (length(missing)) stop(sprintf("%s is missing columns: %s", what, paste(missing, collapse = ", ")))
  build <- attr(x, "genome_build")
  if (is.null(build) || !nzchar(build)) stop(sprintf("%s has no genome_build attribute", what))
  if (anyDuplicated(x$SNP)) {
    stop(sprintf("%s has duplicated SNP ids: %s", what,
                 paste(unique(x$SNP[duplicated(x$SNP)]), collapse = ", ")))
  }
  for (col in c("SNP", "effect_allele", "other_allele", "beta", "se")) {
    if (anyNA(x[[col]])) stop(sprintf("%s has missing values in '%s'", what, col))
  }
  if (any(x$se <= 0)) stop(sprintf("%s has non-positive standard errors", what))
  bad_eaf <- !is.na(x$eaf) & (x$eaf < 0 | x$eaf > 1)
  if (any(bad_eaf)) stop(sprintf("%s has EAF outside [0, 1]", what))
  alleles <- toupper(c(x$effect_allele, x$other_allele))
  if (any(!grepl("^[ACGT]+$", alleles))) stop(sprintf("%s has alleles outside A/C/G/T", what))
  invisible(TRUE)
}

.complement <- function(a) chartr("ACGT", "TGCA", a)

.is_palindrome <- function(a1, a2) {
  nchar(a1) == 1 & nchar(a2) == 1 & .complement(a1) == a2
}

mr_harmonise <- function(exposure, outcome, palindrome_maf_max) {
  if (missing(palindrome_maf_max)) stop("palindrome_maf_max is required (config/mr_methods.yaml)")
  mr_check_sumstats(exposure, "exposure")
  mr_check_sumstats(outcome, "outcome")
  if (!identical(attr(exposure, "genome_build"), attr(outcome, "genome_build"))) {
    stop(sprintf("genome build mismatch: exposure %s vs outcome %s",
                 attr(exposure, "genome_build"), attr(outcome, "genome_build")))
  }

  dropped <- data.frame(SNP = character(), reason = character(), stringsAsFactors = FALSE)
  add_drop <- function(snps, reason) {
    if (length(snps)) dropped <<- rbind(dropped, data.frame(SNP = snps, reason = reason, stringsAsFactors = FALSE))
  }
  add_drop(setdiff(exposure$SNP, outcome$SNP), "absent from outcome")
  add_drop(setdiff(outcome$SNP, exposure$SNP), "absent from exposure")

  m <- merge(exposure[, .mr_required_cols], outcome[, .mr_required_cols],
             by = "SNP", suffixes = c("_exp", "_out"), sort = TRUE)
  e1 <- toupper(m$effect_allele_exp); e2 <- toupper(m$other_allele_exp)
  o1 <- toupper(m$effect_allele_out); o2 <- toupper(m$other_allele_out)
  snv <- nchar(e1) == 1 & nchar(e2) == 1 & nchar(o1) == 1 & nchar(o2) == 1
  pal <- .is_palindrome(e1, e2)

  action <- rep(NA_character_, nrow(m))
  flip_out <- rep(FALSE, nrow(m))  # negate outcome beta and take 1 - EAF

  # Non-palindromic variants: resolved from alleles alone.
  np <- !pal
  direct <- np & e1 == o1 & e2 == o2
  swapped <- np & e1 == o2 & e2 == o1
  sf <- np & snv & !direct & !swapped & e1 == .complement(o1) & e2 == .complement(o2)
  sfs <- np & snv & !direct & !swapped & e1 == .complement(o2) & e2 == .complement(o1)
  action[direct] <- "direct"
  action[swapped] <- "swapped"; flip_out[swapped] <- TRUE
  action[sf] <- "strand_flip"
  action[sfs] <- "strand_flip_swapped"; flip_out[sfs] <- TRUE
  unresolved <- np & is.na(action)
  add_drop(m$SNP[unresolved & snv], "alleles incompatible")
  add_drop(m$SNP[unresolved & !snv], "non-SNV alleles that do not match directly or swapped")

  # Palindromic variants: the allele letters cannot reveal the strand, so the
  # effect-allele frequencies decide, and only when far enough from 0.5.
  same_pair <- pal & ((e1 == o1 & e2 == o2) | (e1 == o2 & e2 == o1))
  add_drop(m$SNP[pal & !same_pair], "alleles incompatible")
  eaf_missing <- same_pair & (is.na(m$eaf_exp) | is.na(m$eaf_out))
  add_drop(m$SNP[eaf_missing], "palindromic, EAF missing")
  cand <- same_pair & !eaf_missing
  maf <- function(f) pmin(f, 1 - f)
  ambiguous <- cand & (maf(m$eaf_exp) > palindrome_maf_max | maf(m$eaf_out) > palindrome_maf_max)
  add_drop(m$SNP[ambiguous], sprintf("palindromic, MAF above %s (strand ambiguous)", palindrome_maf_max))
  ok_pal <- cand & !ambiguous
  # Outcome EAF expressed for the exposure's effect allele as written.
  written_swapped <- ok_pal & e1 == o2
  eaf_out_aligned <- ifelse(written_swapped, 1 - m$eaf_out, m$eaf_out)
  agree <- (m$eaf_exp - 0.5) * (eaf_out_aligned - 0.5) > 0
  # Net flip = (written swapped) XOR (frequencies disagree).
  pal_flip <- ok_pal & (written_swapped != !agree)
  action[ok_pal & !pal_flip] <- "palindrome_direct"
  action[ok_pal & pal_flip] <- "palindrome_flipped"
  flip_out[ok_pal] <- pal_flip[ok_pal]

  keep <- !is.na(action)
  out <- data.frame(
    SNP = m$SNP[keep],
    effect_allele = e1[keep],
    other_allele = e2[keep],
    beta_exp = m$beta_exp[keep],
    se_exp = m$se_exp[keep],
    eaf_exp = m$eaf_exp[keep],
    beta_out = ifelse(flip_out, -m$beta_out, m$beta_out)[keep],
    se_out = m$se_out[keep],
    eaf_out = ifelse(flip_out, 1 - m$eaf_out, m$eaf_out)[keep],
    action = action[keep],
    stringsAsFactors = FALSE
  )
  attr(out, "genome_build") <- attr(exposure, "genome_build")
  n_in <- length(union(exposure$SNP, outcome$SNP))
  if (nrow(out) + nrow(dropped) != n_in) {
    stop(sprintf("harmonisation accounting error: %d kept + %d dropped != %d input variants",
                 nrow(out), nrow(dropped), n_in))
  }
  list(data = out, dropped = dropped,
       counts = c(table(factor(out$action, levels = c("direct", "swapped", "strand_flip",
                                                       "strand_flip_swapped", "palindrome_direct",
                                                       "palindrome_flipped"))),
                  dropped = nrow(dropped)))
}
