# R counterpart of src/common/io.py: analysis inputs must come from data/,
# never from sim/ or tests/, and a missing hand-supplied file stops the run
# with a message naming who provides it. Never creates or fills in a file.

lc_analysis_input <- function(path, provided_by, root = lc_repo_root()) {
  abs <- normalizePath(path, mustWork = FALSE)
  under <- function(dir) startsWith(abs, paste0(normalizePath(dir, mustWork = FALSE), "/"))
  for (d in c("sim", "tests")) {
    if (under(file.path(root, d))) {
      stop(sprintf("%s is under %s/ - simulated/test data must never be read by an analysis pipeline", path, d))
    }
  }
  if (!under(file.path(root, "data"))) stop(sprintf("%s is outside data/ - analysis pipelines read only from data/", path))
  if (!file.exists(abs) || dir.exists(abs)) {
    stop(sprintf(paste("required input %s does not exist. It is provided by: %s.",
                       "Code must not create, scrape or fill in this file - ask for it."),
                 lc_relpath(abs, root), provided_by))
  }
  abs
}
