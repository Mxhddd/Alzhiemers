# LC-First — every figure and table is regenerated from raw inputs by `make all`.
#
#   make setup              install the pinned Python (uv.lock) and R (renv.lock) environments
#   make test               unit tests, Python and R
#   make validate           reproductions of published results (tests/validation/)
#   make provenance-check   Session 0 gate: a run that writes a full provenance header
#   make all                everything that is built so far
#   make lock               re-resolve and re-pin dependencies (then commit the lockfiles)
#   make freeze-prereg NAME=mediation OSF=https://osf.io/xxxxx/
#   make check-prereg NAME=mediation
#
# Module targets are added here in the session that builds each module. Until
# then they fail loudly instead of pretending to run.

SHELL := /bin/bash
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help

PY := uv run --frozen python
RSCRIPT := Rscript

# Modules that are built and gated. Append as each session's gate is passed.
BUILT_MODULES :=

.PHONY: help setup test test-py test-r validate provenance-check all lock \
        freeze-prereg check-prereg clean-results \
        mr baserate mediation dynamics metaanalysis network

help:
	@sed -n '1,14p' Makefile | sed 's/^# \{0,1\}//'

setup:
	uv sync --frozen
	$(RSCRIPT) -e 'renv::restore(prompt = FALSE)'

test: test-py test-r

test-py:
	$(PY) -m pytest tests/unit

test-r:
	$(RSCRIPT) -e 'testthat::test_dir("tests/unit", filter = "provenance", stop_on_failure = TRUE)'

# No validation tests exist until Session 1. pytest exits 5 when it collects
# nothing; that is reported as a failure on purpose: with no validation tests,
# no module is allowed to run.
validate:
	$(PY) -m pytest -m validation tests/validation

provenance-check:
	$(PY) -m src.common.provenance_check --config config/provenance_check.yaml
	$(RSCRIPT) src/common/provenance_check.R config/provenance_check.yaml

all: test provenance-check $(BUILT_MODULES)

lock:
	uv lock
	uv export --frozen --no-hashes --all-groups --format requirements-txt -o requirements.txt
	$(RSCRIPT) -e 'renv::snapshot(prompt = FALSE)'

freeze-prereg:
	@test -n "$(NAME)" -a -n "$(OSF)" || { echo "usage: make freeze-prereg NAME=<module> OSF=https://osf.io/<id>/"; exit 2; }
	$(PY) -m src.common.prereg freeze $(NAME) --document prereg/$(NAME).md --osf-url $(OSF)

check-prereg:
	@test -n "$(NAME)" || { echo "usage: make check-prereg NAME=<module>"; exit 2; }
	$(PY) -m src.common.prereg check $(NAME)

# Removes generated outputs only. Never touches data/raw/.
clean-results:
	find results/figures results/tables results/logs -type f ! -name .gitkeep -delete

mr:           ; @echo "Module 1 (MR) is not built yet — Session 2." >&2; exit 1
baserate:     ; @echo "Module 2 (Markov base-rate) is not built yet — Session 3." >&2; exit 1
metaanalysis: ; @echo "Module 5 (meta-analysis) is not built yet — Session 4." >&2; exit 1
dynamics:     ; @echo "Module 4 (ODE dynamics) is not built yet — Session 5." >&2; exit 1
mediation:    ; @echo "Module 3 (mediation) is not built yet — Session 6." >&2; exit 1
network:      ; @echo "Module 6 (network diffusion) is not built yet — Session 7." >&2; exit 1
