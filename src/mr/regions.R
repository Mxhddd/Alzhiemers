# Gene coordinates for cis windows. Never hardcoded (brief §5): looked up from
# the Ensembl GRCh37 archive via biomaRt, with the Ensembl release logged, and
# cached in data/interim/gene_regions_GRCh37.csv so later runs are offline and
# identical. The cache records the release and retrieval time per gene.

.gene_cache_cols <- c("gene", "chrom", "start", "end", "ensembl_gene_id",
                      "ensembl_release", "retrieved_utc")
.standard_chroms <- c(as.character(1:22), "X", "Y")

mr_gene_cache_path <- function(root = lc_repo_root(), build = "GRCh37") {
  file.path(root, "data", "interim", sprintf("gene_regions_%s.csv", build))
}

mr_read_gene_cache <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, colClasses = c(chrom = "character"))
  miss <- setdiff(.gene_cache_cols, names(x))
  if (length(miss)) stop(sprintf("gene cache %s is missing columns: %s", path, paste(miss, collapse = ", ")))
  if (anyDuplicated(x$gene)) stop(sprintf("gene cache %s has duplicate genes", path))
  x
}

# Query Ensembl GRCh37 for one gene symbol. Network required; errors propagate.
mr_fetch_gene_ensembl <- function(gene) {
  mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", GRCh = 37)
  res <- biomaRt::getBM(
    attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "ensembl_gene_id"),
    filters = "hgnc_symbol", values = gene, mart = mart)
  res <- res[res$chromosome_name %in% .standard_chroms, , drop = FALSE]
  if (nrow(res) != 1) {
    stop(sprintf("Ensembl GRCh37 returned %d records on standard chromosomes for %s; resolve by hand",
                 nrow(res), gene))
  }
  marts <- biomaRt::listMarts(host = "https://grch37.ensembl.org")
  release <- marts$version[marts$biomart == "ENSEMBL_MART_ENSEMBL"]
  if (length(release) != 1) stop("could not determine the Ensembl GRCh37 release")
  data.frame(gene = gene, chrom = res$chromosome_name, start = res$start_position,
             end = res$end_position, ensembl_gene_id = res$ensembl_gene_id,
             ensembl_release = release, retrieved_utc = lc_iso_utc(lc_utc_now()),
             stringsAsFactors = FALSE)
}

# Region = gene body +/- window_kb. Uses the cache; on a miss, queries Ensembl
# (unless offline = TRUE, which stops instead) and appends to the cache.
mr_gene_region <- function(gene, window_kb, run = NULL, root = lc_repo_root(), build = "GRCh37",
                           offline = FALSE, fetch = mr_fetch_gene_ensembl) {
  if (build != "GRCh37") stop("gene_region() only supports GRCh37 (the outcome build)")
  if (!is.numeric(window_kb) || length(window_kb) != 1 || window_kb < 0) stop("window_kb must be a non-negative number")
  path <- mr_gene_cache_path(root, build)
  cache <- mr_read_gene_cache(path)
  hit <- if (is.null(cache)) NULL else cache[cache$gene == gene, , drop = FALSE]
  if (is.null(hit) || nrow(hit) == 0) {
    if (offline) stop(sprintf("gene %s is not in %s and offline = TRUE", gene, path))
    hit <- fetch(gene)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    cache <- rbind(cache, hit[, .gene_cache_cols])
    utils::write.csv(cache, path, row.names = FALSE)
  }
  if (!is.null(run)) {
    lc_log(run, sprintf("gene_region %s: chr%s:%d-%d (%s, Ensembl GRCh37 release %s, retrieved %s)",
                        gene, hit$chrom, hit$start, hit$end, hit$ensembl_gene_id,
                        hit$ensembl_release, hit$retrieved_utc))
  }
  list(gene = gene, chrom = as.character(hit$chrom),
       start = max(1, hit$start - window_kb * 1000), end = hit$end + window_kb * 1000,
       gene_start = hit$start, gene_end = hit$end, build = build,
       ensembl_release = hit$ensembl_release)
}
