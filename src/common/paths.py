"""Canonical repository paths. Every module resolves files through these so
that the input guards in ``io.py`` can tell real data from simulated data."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

CONFIG_DIR = REPO_ROOT / "config"
PREREG_DIR = REPO_ROOT / "prereg"
DATA_DIR = REPO_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
SIM_DIR = REPO_ROOT / "sim"
TESTS_DIR = REPO_ROOT / "tests"
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"
LOGS_DIR = RESULTS_DIR / "logs"

# Lockfiles whose hashes go into every run header: together they pin the
# whole software environment, not just the handful of packages a module names.
LOCKFILES = (REPO_ROOT / "uv.lock", REPO_ROOT / "renv.lock")
