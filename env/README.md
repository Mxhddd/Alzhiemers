# env/

Environment notes. The pinned environments themselves live at the repo root:
`uv.lock` / `requirements.txt` (Python) and `renv.lock` (R). The SHA-256 of both
lockfiles is written into every run log header.

## Versions used to build Session 0

| Tool | Version | Notes |
|---|---|---|
| Python | 3.11.x (`.python-version`) | `requires-python = ">=3.11,<3.12"` — one interpreter for reproducibility |
| uv | 0.11 | `uv sync --frozen` installs exactly `uv.lock` |
| R | 4.3.3 | Ubuntu 24.04 `r-base-core` |
| renv | 1.0.3 | |

## Known limitation of the build container

The cloud container used for Session 0 could not reach CRAN (network policy),
so R packages were installed from Ubuntu's `r-cran-*` builds, which are the
CRAN releases listed in `renv.lock`. On a normal machine `renv::restore()`
fetches the same versions from CRAN. Later sessions need Bioconductor/CRAN
packages (TwoSampleMR, coloc, metafor, ...) and therefore network access to
CRAN — or the matching apt packages where they exist.
