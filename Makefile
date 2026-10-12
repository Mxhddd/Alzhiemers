# LC-First — every figure and table is regenerated from raw inputs by `make all`.
#
#   make setup              install the pinned Python (uv.lock) and R (renv.lock) environments
#   make test               unit tests, Python and R
#   make validate           reproductions of published results (tests/validation/)
#   make check-validated MODULE=mr   is the module's last validation still valid?
#   make provenance-check   Session 0 gate: a run that writes a full provenance header
#   make all                everything that is built so far
#   make lock               re-resolve and re-pin dependencies (then commit the lockfiles)
#   make freeze-prereg NAME=mediation OSF=https://osf.io/xxxxx/
#   make check-prereg NAME=mediation
#
# Module targets are added here in the session that builds each module. Until
# then they fail loudly instead of pretending to run.

SHELL := /bin/bash
# UTF-8 everywhere: in a C locale R misreads non-ASCII text.
export LC_ALL := C.UTF-8
.SHELLFLAGS := -eu -o pipefail -c
.DEFAULT_GOAL := help

PY := uv run --frozen python
RSCRIPT := Rscript

# Modules that are built and gated. Append as each session's gate is passed.
BUILT_MODULES :=

.PHONY: help setup test test-py test-r validate validate-mr check-validated provenance-check all lock \
        freeze-prereg check-prereg clean-results \
        mr baserate mediation dynamics metaanalysis network

help:
	@sed -n '1,15p' Makefile | sed 's/^# \{0,1\}//'

setup:
	uv sync --frozen
	$(RSCRIPT) -e 'renv::restore(prompt = FALSE)'

test: test-py test-r

test-py:
	$(PY) -m pytest tests/unit

test-r:
	$(RSCRIPT) -e 'testthat::test_dir("tests/unit", stop_on_failure = TRUE, stop_on_warning = TRUE)'

# Reproductions of published results. Each module's run records PASSED/FAILED
# plus a code fingerprint in results/validation/<module>.status.yaml; analysis
# entry points refuse to run unless it passed on the current code.
# A target whose published numbers have not been transcribed FAILS (never skips).
validate: validate-mr

validate-mr:
	$(RSCRIPT) tests/validation/run_validation.R mr

check-validated:
	@test -n "$(MODULE)" || { echo "usage: make check-validated MODULE=<module>"; exit 2; }
	$(PY) -m src.common.validation check $(MODULE)

provenance-check:
	$(PY) -m src.common.provenance_check --config config/provenance_check.yaml
	$(RSCRIPT) src/common/provenance_check.R config/provenance_check.yaml

all: test validate provenance-check $(BUILT_MODULES)

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
	find results/figures results/tables results/logs results/validation -type f ! -name .gitkeep -delete

mr:           ; @echo "Module 1 (MR) is not built yet — Session 2." >&2; exit 1
baserate:     ; @echo "Module 2 (Markov base-rate) is not built yet — Session 3." >&2; exit 1
metaanalysis: ; @echo "Module 5 (meta-analysis) is not built yet — Session 4." >&2; exit 1
dynamics:     ; @echo "Module 4 (ODE dynamics) is not built yet — Session 5." >&2; exit 1
mediation:    ; @echo "Module 3 (mediation) is not built yet — Session 6." >&2; exit 1
network:      ; @echo "Module 6 (network diffusion) is not built yet — Session 7." >&2; exit 1
