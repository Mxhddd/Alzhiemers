# LC-First — every figure and table is regenerated from raw inputs by `make all`.
#
#   make setup              install the pinned Python (uv.lock) and R (renv.lock) environments
#   make test               unit tests, Python and R
#   make validate           reproductions of published results (tests/validation/)
#   make check-validated MODULE=mr   is the module's last validation still valid?
#   make provenance-check   Session 0 gate: a run that writes a full provenance header
#   make all                everything that is built so far
#   make todo               list every input the researcher still has to supply
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

# Every module, in dependency order (metaanalysis triangulates with mr). Each
# refuses to run until its validation passed and its prereg is frozen.
BUILT_MODULES := mr baserate metaanalysis dynamics mediation network

.PHONY: help todo setup test test-py test-r test-slow validate validate-mr validate-baserate validate-metaanalysis validate-dynamics validate-mediation validate-network check-validated provenance-check all lock \
        freeze-prereg check-prereg clean-results \
        mr baserate mediation dynamics metaanalysis network

help:
	@sed -n '1,16p' Makefile | sed 's/^# \{0,1\}//'

setup:
	uv sync --frozen
	$(RSCRIPT) -e 'renv::restore(prompt = FALSE)'

test: test-py test-r

test-py:
	$(PY) -m pytest -m "not slow" tests/unit

# Minutes-long tests (MCMC). Never skipped: `make all` runs them.
test-slow:
	$(PY) -m pytest -m slow tests/unit

test-r:
	$(RSCRIPT) -e 'testthat::test_dir("tests/unit", stop_on_failure = TRUE, stop_on_warning = TRUE)'

# Reproductions of published results. Each module's run records PASSED/FAILED
# plus a code fingerprint in results/validation/<module>.status.yaml; analysis
# entry points refuse to run unless it passed on the current code.
# A target whose published numbers have not been transcribed FAILS (never skips).
validate: validate-mr validate-baserate validate-metaanalysis validate-dynamics validate-mediation validate-network

validate-mr:
	$(RSCRIPT) tests/validation/run_validation.R mr

validate-metaanalysis:
	$(RSCRIPT) tests/validation/run_validation.R metaanalysis

validate-dynamics:
	$(PY) -m tests.validation.run_validation_py dynamics

validate-mediation:
	$(PY) -m tests.validation.run_validation_py mediation

validate-network:
	$(PY) -m tests.validation.run_validation_py network

validate-baserate:
	$(PY) -m tests.validation.run_validation_py baserate

todo:
	@$(PY) -m src.common.todo

check-validated:
	@test -n "$(MODULE)" || { echo "usage: make check-validated MODULE=<module>"; exit 2; }
	$(PY) -m src.common.validation check $(MODULE)

provenance-check:
	$(PY) -m src.common.provenance_check --config config/provenance_check.yaml
	$(RSCRIPT) src/common/provenance_check.R config/provenance_check.yaml

all: test test-slow validate provenance-check $(BUILT_MODULES)

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

# Module 1. Refuses to run unless `make validate` passed on the current code
# and prereg/mr.md is frozen; then draws the forest plot.
mr:
	$(RSCRIPT) src/mr/run_mr.R
	$(PY) -m src.mr.forest_plot --summary results/tables/mr_summary.csv
# Module 2. Refuses to run unless validated on the current code and
# prereg/baserate.md is frozen; needs the researcher's hand-transcribed files.
baserate:
	$(PY) -m src.baserate.run_baserate
# Module 5. Gated on validation + frozen prereg/metaanalysis.md; the
# anticholinergic positive control must be recovered before anything else.
metaanalysis:
	$(RSCRIPT) src/metaanalysis/run_meta.R
	$(PY) -m src.metaanalysis.plots
# Module 4. Gated on validation (incl. the researcher's analytic fixed point)
# + frozen prereg/dynamics.md; the LC model equations are the researcher's.
dynamics:
	$(PY) -m src.dynamics.run_dynamics
# Module 3 (Killer 2). Real data only; gated on validation + frozen DAG prereg.
mediation:
	$(PY) -m src.mediation.run_mediation
# Module 6. Real data only; gated on validation (spin test vs known null) + frozen prereg.
network:
	$(PY) -m src.network.run_network
