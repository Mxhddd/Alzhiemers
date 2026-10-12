# Instrument selection: p-value, MAF and F filters, then LD clumping.
# Each filter's drops are returned with reasons (brief: "F > 10, report drops").

# Greedy clumping (as PLINK --clump): repeatedly take the most significant
# remaining variant, remove every other remaining variant within window_kb
# whose r^2 with it exceeds r2_max. `r2_fun(snps)` returns a symmetric r^2
# matrix with dimnames = the SNPs present in the LD reference; variants absent
# from the reference are dropped (their LD is unknown, so independence cannot
# be checked).
mr_clump <- function(x, r2_fun, r2_max, window_kb) {
  dropped <- data.frame(SNP = character(), reason = character(), stringsAsFactors = FALSE)
  if (nrow(x) <= 1) return(list(data = x, dropped = dropped))
  r2 <- r2_fun(x$SNP)
  if (!is.matrix(r2) || !identical(rownames(r2), colnames(r2))) stop("r2_fun must return a square matrix with matching dimnames")
  absent <- setdiff(x$SNP, rownames(r2))
  if (length(absent)) dropped <- rbind(dropped, data.frame(SNP = absent, reason = "absent from LD reference panel"))
  x <- x[x$SNP %in% rownames(r2), , drop = FALSE]
  x <- x[order(x$p, -abs(x$beta / x$se)), , drop = FALSE]
  remaining <- x$SNP
  kept <- character()
  while (length(remaining)) {
    lead <- remaining[1]
    kept <- c(kept, lead)
    others <- remaining[-1]
    pos_lead <- x$pos[x$SNP == lead]
    near <- abs(x$pos[match(others, x$SNP)] - pos_lead) <= window_kb * 1000
    linked <- others[near & r2[lead, others] > r2_max]
    if (length(linked)) {
      dropped <- rbind(dropped, data.frame(SNP = linked, reason = sprintf("LD r2 > %s with %s", r2_max, lead)))
    }
    remaining <- setdiff(others, linked)
  }
  list(data = x[x$SNP %in% kept, , drop = FALSE], dropped = dropped)
}

# r^2 from a PLINK 1.9 binary reference panel. Requires the `plink` binary.
mr_plink_r2_fun <- function(bfile, plink = Sys.which("plink")) {
  if (!nzchar(plink)) stop("PLINK 1.9 binary not found on PATH; LD clumping cannot run")
  for (ext in c(".bed", ".bim", ".fam")) {
    if (!file.exists(paste0(bfile, ext))) stop(sprintf("LD reference %s%s not found", bfile, ext))
  }
  function(snps) {
    tmp <- tempfile()
    writeLines(snps, paste0(tmp, ".snps"))
    out <- suppressWarnings(system2(plink, c("--bfile", shQuote(bfile), "--extract", paste0(tmp, ".snps"),
                                             "--r2", "square", "--write-snplist", "--out", tmp),
                                    stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(out, "status"))) stop(paste(c("plink failed:", out), collapse = "\n"))
    ids <- readLines(paste0(tmp, ".snplist"))
    m <- as.matrix(utils::read.table(paste0(tmp, ".ld")))
    dimnames(m) <- list(ids, ids)
    m
  }
}

mr_build_instruments <- function(x, instr, r2_fun) {
  dropped <- data.frame(SNP = character(), reason = character(), stringsAsFactors = FALSE)
  step <- function(keep, reason) {
    if (any(!keep)) dropped <<- rbind(dropped, data.frame(SNP = x$SNP[!keep], reason = reason, stringsAsFactors = FALSE))
    x <<- x[keep, , drop = FALSE]
  }
  step(x$p <= instr$p_max, sprintf("p > %g", instr$p_max))
  step(pmin(x$eaf, 1 - x$eaf) >= instr$maf_min, sprintf("MAF < %g", instr$maf_min))
  x$f_stat <- mr_f_stat(x$beta, x$se)
  step(x$f_stat > instr$f_stat_min, sprintf("F <= %g (weak instrument)", instr$f_stat_min))
  cl <- mr_clump(x, r2_fun, instr$ld_r2_max, instr$ld_window_kb)
  dropped <- rbind(dropped, cl$dropped)
  list(data = cl$data, dropped = dropped)
}
