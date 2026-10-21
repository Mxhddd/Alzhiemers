"""Download the public data for Module 3's two published-reproduction checks.

    python -m src.mediation.fetch_validation_data

Both are the authors' own replication archives on Harvard Dataverse, released
under CC0 1.0, used ONLY to validate the estimator (never as Module 3 data):

* Imai, Keele & Tingley (2010) Psychol Methods 15:309, JOBS II, Table 1.
  Dataverse doi:10.7910/DVN/UMEYXD, file Jobs-NoMiss-Cont (Stata .dta).
* Imai, Keele & Yamamoto (2010) Stat Sci 25:51, Nelson et al. framing data,
  Table 3 and Section 6.2. Dataverse doi:10.7910/DVN/APNMP7, nelson_data.RData.

Each file is checked against a pinned SHA-256 before use, then converted to
CSV with base R (`foreign::read.dta` / `load`), keeping the authors' factor
labels. Output: data/raw/validation/jobs_ii_ikt2010.csv and
data/raw/validation/nelson_iky2010.csv, plus MANIFEST_mediation.tsv.
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from src.common.paths import REPO_ROOT

OUT = REPO_ROOT / "data" / "raw" / "validation"
FILES = [
    {"out": "jobs_ii_ikt2010.csv",
     "url": "https://dataverse.harvard.edu/api/access/datafile/1576667?format=original",
     "sha256": "c152ee529c51ff6303cdd562238dd1a89603fbb5e7b7125767055e3e58fe21df",
     "r": 'library(foreign); d <- read.dta(commandArgs(TRUE)[1], convert.underscore = FALSE); '
          'write.csv(d, commandArgs(TRUE)[2], row.names = FALSE)',
     "source": "Imai, Keele & Tingley 2010, doi:10.7910/DVN/UMEYXD (Jobs-NoMiss-Cont.dta), CC0"},
    {"out": "nelson_iky2010.csv",
     "url": "https://dataverse.harvard.edu/api/access/datafile/1545812",
     "sha256": "6b566febe6620553943ed839ab40bc67874bc9b7cb46372dc42e395ccb9c4d69",
     "r": 'e <- new.env(); load(commandArgs(TRUE)[1], envir = e); '
          'write.csv(e$kkk, commandArgs(TRUE)[2], row.names = FALSE)',
     "source": "Imai, Keele & Yamamoto 2010, doi:10.7910/DVN/APNMP7 (nelson_data.RData), CC0"},
]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    with tempfile.TemporaryDirectory() as td:
        for f in FILES:
            raw = Path(td) / (f["out"] + ".download")
            # curl, not urllib: it honours the system proxy / CA settings the same way everywhere
            subprocess.run(["curl", "-sS", "-f", "-L", "-m", "300", "-o", str(raw), f["url"]], check=True)
            digest = hashlib.sha256(raw.read_bytes()).hexdigest()
            if digest != f["sha256"]:
                raise RuntimeError(f"{f['url']}: sha256 {digest} != pinned {f['sha256']}; not used")
            dest = OUT / f["out"]
            subprocess.run(["Rscript", "--vanilla", "-e", f["r"], str(raw), str(dest)], check=True,
                           cwd=td, capture_output=True, text=True)
            rows.append("\t".join([str(dest.relative_to(REPO_ROOT)), f["url"], f["sha256"],
                                   hashlib.sha256(dest.read_bytes()).hexdigest(),
                                   datetime.now(timezone.utc).isoformat(timespec="seconds"), f["source"]]))
            print(f"wrote {dest.relative_to(REPO_ROOT)}")
    (OUT / "MANIFEST_mediation.tsv").write_text(
        "path\turl\tsha256_download\tsha256_csv\tutc\tsource\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
