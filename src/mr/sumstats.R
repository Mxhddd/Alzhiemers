# Reading GWAS / QTL summary statistics into the harmonise() schema.
#
# Every source has a column map in config/mr_targets.yaml (our name -> file
# header). Required: SNP, chrom, pos, effect_allele, other_allele, eaf, and
# either beta + se, or z + n (eQTLGen gives only Z-scores). Optional: p,
# gene (QTL files: rows for many genes), n.
#
# Rows are never silently removed: anything dropped is returned in `dropped`
# with a reason, for the caller to log.

.ss_required <- c("SNP", "chrom", "pos", "effect_allele", "other_allele", "eaf")

# eQTLGen-style Z-score -> per-SD beta and SE (Zhu et al. 2016, Nat Genet,
# SMR methods): beta = z / sqrt(2p(1-p)(n + z^2)), se = 1 / sqrt(2p(1-p)(n + z^2)).
mr_beta_from_z <- function(z, eaf, n) {
  d <- sqrt(2 * eaf * (1 - eaf) * (n + z^2))
  list(beta = z / d, se = 1 / d)
}

mr_read_sumstats <- function(path, columns, build, provided_by, region = NULL, gene = NULL,
                             root = lc_repo_root()) {
  if (is.null(path)) stop(sprintf("no file configured (provided by: %s)", provided_by))
  if (is.null(columns)) stop(sprintf("no column map configured for %s", path))
  if (is.null(build)) stop(sprintf("no genome build configured for %s", path))
  abs <- lc_analysis_input(file.path(root, path), provided_by, root = root)
  has_beta <- all(c("beta", "se") %in% names(columns))
  has_z <- all(c("z", "n") %in% names(columns))
  if (!has_beta && !has_z) stop(sprintf("column map for %s needs beta+se or z+n", path))
  miss <- setdiff(.ss_required, names(columns))
  if (length(miss)) stop(sprintf("column map for %s lacks: %s", path, paste(miss, collapse = ", ")))

  # data.table::fread reads .gz files through R.utils (in renv.lock since 2026-10-21).
  if (grepl("\\.(gz|bz2)$", abs) && !requireNamespace("R.utils", quietly = TRUE)) {
    stop(sprintf("%s is compressed and the R.utils package is not installed (renv::restore())", path))
  }
  header <- names(data.table::fread(abs, nrows = 0))
  absent <- setdiff(unlist(columns), header)
  if (length(absent)) stop(sprintf("%s has no column(s) %s (header: %s)", path,
                                   paste(absent, collapse = ", "), paste(header, collapse = ", ")))
  x <- data.table::fread(abs, select = unname(unlist(columns)), data.table = FALSE,
                         colClasses = list(character = unname(unlist(columns[c("SNP", "chrom")]))))
  names(x) <- names(columns)[match(names(x), unlist(columns))]
  n_read <- nrow(x)
  dropped <- data.frame(SNP = character(), reason = character(), stringsAsFactors = FALSE)
  drop <- function(keep, reason) {
    if (any(!keep)) dropped <<- rbind(dropped, data.frame(SNP = x$SNP[!keep], reason = reason, stringsAsFactors = FALSE))
    x <<- x[keep, , drop = FALSE]
  }

  if (!is.null(gene)) {
    if (!"gene" %in% names(x)) stop(sprintf("%s: a gene filter was requested but no 'gene' column is mapped", path))
    x <- x[x$gene == gene, , drop = FALSE]   # rows for other genes are not drops
  }
  x$chrom <- sub("^chr", "", x$chrom, ignore.case = TRUE)
  if (!is.null(region)) {
    x <- x[x$chrom == region$chrom & x$pos >= region$start & x$pos <= region$end, , drop = FALSE]
  }
  x$effect_allele <- toupper(x$effect_allele); x$other_allele <- toupper(x$other_allele)
  drop(!is.na(x$eaf), "EAF missing")
  if (!has_beta) {
    drop(!is.na(x$z) & !is.na(x$n), "Z or N missing")
    d <- mr_beta_from_z(x$z, x$eaf, x$n)
    x$beta <- d$beta; x$se <- d$se
  }
  drop(!is.na(x$beta) & !is.na(x$se) & x$se > 0, "beta/SE missing or SE <= 0")
  if (!"p" %in% names(x)) x$p <- 2 * stats::pnorm(-abs(x$beta / x$se))
  drop(!duplicated(x$SNP) & !duplicated(x$SNP, fromLast = TRUE), "duplicated SNP id (multi-allelic or repeated row)")
  attr(x, "genome_build") <- build
  list(data = x, dropped = dropped, n_read = n_read)
}
