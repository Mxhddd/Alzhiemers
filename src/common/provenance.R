# Provenance and run logging for the R modules (Module 1 MR, Module 5 meta-analysis).
#
# Writes the same log format as src/common/runlog.py: a YAML header between
# `---` lines, log lines, and a `--- footer` block. Usage:
#
#   source("src/common/provenance.R")
#   lc_with_run("mr", "config/mr_targets.yaml", seed = cfg$seed,
#               packages = c("TwoSampleMR", "coloc"), function(run) {
#     lc_log(run, "fetching exposures")
#     lc_drop(run, stage = "instruments", item = "rs123", reason = "F < 10")
#     lc_record_output(run, "results/tables/mr_summary.csv")
#   })
#
# Errors are never caught: if the body fails, the footer says FAILED and the
# error propagates. The pre-registration gate is NOT reimplemented here; R
# calls the single Python implementation so the two languages cannot diverge.

lc_repo_root <- function(start = getwd()) {
  dir <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(dir, "PROJECT_BRIEF.md")) && dir.exists(file.path(dir, ".git"))) {
      return(dir)
    }
    parent <- dirname(dir)
    if (identical(parent, dir)) stop("not inside the lc-first repository (no PROJECT_BRIEF.md + .git found)")
    dir <- parent
  }
}

lc_relpath <- function(path, root) {
  prefix <- paste0(root, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1) else path
}

lc_git <- function(args, root) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(root), args), stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(sprintf("git %s failed: %s", paste(args, collapse = " "), paste(out, collapse = "\n")))
  }
  out
}

lc_git_state <- function(root = lc_repo_root()) {
  commit <- lc_git(c("rev-parse", "HEAD"), root)
  status <- lc_git(c("status", "--porcelain=v1", "--untracked-files=all"), root)
  status <- status[nzchar(trimws(status))]
  files <- substring(status, 4)
  list(commit = commit[[1]], dirty = length(files) > 0, uncommitted_files = as.list(files))
}

lc_sha256 <- function(path) {
  if (!file.exists(path) || dir.exists(path)) stop(sprintf("cannot hash %s: not a file", path))
  digest::digest(file = path, algo = "sha256")
}

lc_utc_now <- function() {
  t <- Sys.time()
  attr(t, "tzone") <- "UTC"
  t
}

lc_iso_utc <- function(t) format(t, "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")

lc_package_versions <- function(packages) {
  versions <- list(R = paste(R.version$major, R.version$minor, sep = "."))
  for (p in packages) {
    # packageVersion() errors if the package is absent: we want that.
    versions[[p]] <- as.character(utils::packageVersion(p))
  }
  versions
}

lc_start_run <- function(module, config_path, seed, packages, preregistration = NULL,
                         root = lc_repo_root(), logs_dir = file.path(root, "results", "logs")) {
  if (missing(seed)) stop("seed is required: give the config's seed, or NULL for a deterministic run")
  if (!nzchar(module) || grepl("/", module, fixed = TRUE)) stop(sprintf("bad module name: %s", module))
  started <- lc_utc_now()
  git <- lc_git_state(root)
  config_abs <- normalizePath(config_path, mustWork = TRUE)
  config_rel <- lc_relpath(config_abs, root)
  lockfiles <- list()
  for (lf in c("uv.lock", "renv.lock")) lockfiles[[lf]] <- lc_sha256(file.path(root, lf))

  header <- list(
    lc_first_run_header = 1L,
    module = module,
    language = "R",
    started_utc = lc_iso_utc(started),
    git = git,
    config = list(path = config_rel, sha256 = lc_sha256(config_abs)),
    lockfiles = lockfiles,
    seed = seed,
    packages = lc_package_versions(packages),
    platform = R.version$platform,
    command = as.list(commandArgs()),
    preregistration = preregistration
  )

  log_dir <- file.path(logs_dir, module)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  stamp <- format(started, "%Y%m%dT%H%M%OS6Z", tz = "UTC")
  stamp <- sub(".", "", stamp, fixed = TRUE)
  log_path <- file.path(log_dir, sprintf("%s_%s.log", stamp, substr(git$commit, 1, 7)))
  if (file.exists(log_path)) stop(sprintf("log %s already exists; refusing to overwrite", log_path))
  writeLines(c("---", sub("\n$", "", yaml::as.yaml(header, handlers = list(
    logical = function(x) { structure(ifelse(x, "true", "false"), class = "verbatim") }
  ))), "---"), log_path)

  run <- new.env()
  run$module <- module
  run$log_path <- log_path
  run$header <- header
  run$root <- root
  run$drops <- data.frame(stage = character(), item = character(), reason = character(),
                          stringsAsFactors = FALSE)
  run$outputs <- list()
  run$closed <- FALSE
  lc_log(run, paste("log:", log_path))
  if (isTRUE(git$dirty)) {
    lc_log(run, paste0("WORKING TREE IS DIRTY - this result is not reproducible from commit ",
                       substr(git$commit, 1, 7), " alone. Uncommitted: ",
                       paste(unlist(git$uncommitted_files), collapse = ", ")),
           level = "WARNING")
  }
  run
}

lc_log <- function(run, msg, level = "INFO") {
  line <- sprintf("%s %s %s", format(lc_utc_now(), "%Y-%m-%d %H:%M:%S", tz = "UTC"), level, msg)
  cat(line, "\n", sep = "", file = run$log_path, append = TRUE)
  message(line)
  invisible(NULL)
}

lc_drop <- function(run, stage, item, reason) {
  for (v in list(stage, item, reason)) {
    if (length(v) != 1 || is.na(v) || !nzchar(trimws(v))) stop("lc_drop needs non-empty stage, item and reason")
  }
  run$drops[nrow(run$drops) + 1, ] <- c(stage, item, reason)
  lc_log(run, sprintf("DROPPED [%s] %s - %s", stage, item, reason), level = "WARNING")
}

lc_record_output <- function(run, path) {
  abs <- normalizePath(path, mustWork = TRUE)
  key <- lc_relpath(abs, run$root)
  run$outputs[[key]] <- lc_sha256(abs)
  lc_log(run, sprintf("OUTPUT %s sha256=%s", key, run$outputs[[key]]))
}

lc_empty_map <- function() structure(list(), names = character())  # YAML {} like Python

lc_drop_counts <- function(run) {
  out <- lc_empty_map()
  for (st in unique(run$drops$stage)) {
    rows <- run$drops[run$drops$stage == st, , drop = FALSE]
    tab <- table(rows$reason)
    out[[st]] <- as.list(stats::setNames(as.integer(tab), names(tab)))
  }
  out
}

lc_write_drop_report <- function(run, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c(sprintf("# Dropped items - %s", run$module), "",
             sprintf("Run log: `%s`  ", lc_relpath(normalizePath(run$log_path), run$root)),
             sprintf("Commit: `%s`  ", run$header$git$commit),
             sprintf("Total dropped: **%d**", nrow(run$drops)), "")
  if (nrow(run$drops) > 0) {
    lines <- c(lines, "| stage | item | reason |", "|---|---|---|",
               sprintf("| %s | %s | %s |", run$drops$stage, run$drops$item,
                       gsub("|", "/", run$drops$reason, fixed = TRUE)))
  } else {
    lines <- c(lines, "Nothing was dropped.")
  }
  writeLines(lines, path)
  lc_record_output(run, path)
  invisible(path)
}

lc_finish_run <- function(run, status, error = NULL) {
  if (isTRUE(run$closed)) stop("run already closed")
  run$closed <- TRUE
  footer <- list(
    finished_utc = lc_iso_utc(lc_utc_now()),
    status = status,
    error = error,
    n_dropped = nrow(run$drops),
    dropped_by_stage_and_reason = lc_drop_counts(run),
    outputs_sha256 = if (length(run$outputs)) run$outputs else lc_empty_map()
  )
  cat("--- footer\n", yaml::as.yaml(footer), "---\n", sep = "", file = run$log_path, append = TRUE)
  invisible(run)
}

# Run `body(run)` inside a logged run. The footer is always written; errors propagate.
lc_with_run <- function(module, config_path, seed, packages, body, ...) {
  run <- lc_start_run(module, config_path, seed, packages, ...)
  err <- NULL
  result <- withCallingHandlers(
    body(run),
    error = function(e) {
      err <<- conditionMessage(e)
      lc_log(run, paste("RUN FAILED:", err), level = "ERROR")
      lc_finish_run(run, "FAILED", err)
    }
  )
  lc_finish_run(run, "completed")
  invisible(run)
}

# Pre-registration gate: delegates to the Python implementation (single source of truth).
lc_require_preregistration <- function(name, root = lc_repo_root()) {
  out <- suppressWarnings(system2("uv", c("--directory", shQuote(root), "run", "--frozen", "python", "-m",
                                                 "src.common.prereg", "check", name),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop(paste(c(sprintf("pre-registration gate for '%s' FAILED:", name), out), collapse = "\n"))
  }
  yaml::read_yaml(file.path(root, "prereg", paste0(name, ".freeze.yaml")))
}

lc_read_header <- function(log_path) {
  lines <- readLines(log_path)
  if (length(lines) == 0 || lines[[1]] != "---") stop(sprintf("%s: no provenance header", log_path))
  end <- which(lines == "---")[-1][1]
  if (is.na(end)) stop(sprintf("%s: unterminated provenance header", log_path))
  yaml::yaml.load(paste(lines[2:(end - 1)], collapse = "\n"))
}
