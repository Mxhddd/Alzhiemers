"""Package the pre-registrations for OSF upload (never freezes anything).

    uv run python -m src.common.prereg_package      (or: make prereg-package)

Writes results/osf_package/<UTC date>_<commit>/:
  <name>.md                     each pre-registration document, as committed
  config/...                    the configuration files each document fixes
  MANIFEST.tsv                  name, status, SHA-256 of every file, commit
  UPLOAD_INSTRUCTIONS.md        what the researcher does next
and a .zip of the same folder. Refuses to package from a dirty tree, so
the hashes match the commit. Freezing stays the researcher's step:
`make freeze-prereg NAME=<name> OSF=<url>` after uploading the document.
"""

from __future__ import annotations

import shutil
import sys
import time
from pathlib import Path

from src.common.paths import PREREG_DIR, REPO_ROOT, RESULTS_DIR
from src.common.prereg import BLOCK_MARKER
from src.common.provenance import git_state, sha256_file

CONFIGS = {
    "mr": ["config/mr_targets.yaml", "config/mr_methods.yaml"],
    "baserate": ["config/baserate.yaml"],
    "metaanalysis": ["config/metaanalysis.yaml"],
    "dynamics": ["config/dynamics.yaml", "config/dynamics_params.yaml"],
    "mediation": ["config/mediation.yaml"],
    "network": ["config/network.yaml"],
}


def main() -> int:
    git = git_state(REPO_ROOT)
    if git.dirty:
        raise SystemExit("working tree is dirty: commit first so the package matches a commit")
    stamp = time.strftime("%Y%m%d", time.gmtime())
    out = RESULTS_DIR / "osf_package" / f"{stamp}_{git.commit[:7]}"
    if out.exists():
        shutil.rmtree(out)
    (out / "config").mkdir(parents=True)
    rows = ["name\tstatus\tfile\tsha256"]
    ready, blocked = [], []
    for name, cfgs in CONFIGS.items():
        doc = PREREG_DIR / f"{name}.md"
        if not doc.exists():
            rows.append(f"{name}\tNOT WRITTEN\t-\t-")
            continue
        frozen = (PREREG_DIR / f"{name}.freeze.yaml").exists()
        status = ("FROZEN" if frozen else "BLOCKED (DO NOT FREEZE marker)" if BLOCK_MARKER in doc.read_text()
                  else "READY TO FREEZE")
        (ready if status == "READY TO FREEZE" else blocked).append((name, status))
        shutil.copy2(doc, out / doc.name)
        rows.append(f"{name}\t{status}\t{doc.name}\t{sha256_file(doc)}")
        for c in cfgs:
            src = REPO_ROOT / c
            shutil.copy2(src, out / "config" / src.name)
            rows.append(f"{name}\t(config)\t{c}\t{sha256_file(src)}")
    (out / "MANIFEST.tsv").write_text("\n".join(rows) + f"\n# commit {git.commit}\n")
    lines = [f"# OSF upload package — commit {git.commit}", "",
             "Built by `make prereg-package`. Nothing here has been frozen.",
             "Freezing is the researcher's step. The agent that wrote the criteria must not lock them.", "",
             "## Ready to freeze", ""]
    lines += [f"- `{n}.md`" for n, _ in ready] or ["- (none)"]
    lines += ["", "## Not ready", ""] + ([f"- `{n}.md`: {s}" for n, s in blocked] or ["- (none)"])
    lines += ["", "## Steps for each ready document", "",
              "1. Read it, including its change log and the config files in `config/` that it fixes.",
              "2. On OSF, create a registration (or upload the file to an OSF project) and note its URL.",
              "3. In the repository, at the SAME commit, run "
              "`make freeze-prereg NAME=<name> OSF=https://osf.io/<id>/`, then commit "
              "`prereg/<name>.freeze.yaml`.",
              "4. `make check-prereg NAME=<name>` must pass. After that, the module may run on real data.",
              "", "SHA-256 values for every file are in `MANIFEST.tsv`. "
              "The freeze tool records the document's hash, and any later edit breaks the gate."]
    (out / "UPLOAD_INSTRUCTIONS.md").write_text("\n".join(lines) + "\n")
    z = shutil.make_archive(str(out), "zip", root_dir=out.parent, base_dir=out.name)
    print(f"{out.relative_to(REPO_ROOT)} (+ {Path(z).name}): ready {[n for n, _ in ready]}, not ready {blocked}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
