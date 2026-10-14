# Stage 4 — researcher's update of 2026-10-13 (network, Braak 2011, criterion 2, Modules 5 and 1)

Done 2026-10-14.

- **GitHub:** the push works. Expected head f1fea88 confirmed on the remote. Every later commit is pushed.
- **Network audit:** `manuscript/notes/network_audit.md`.
  - Nothing is blocked by the environment allowlist.
  - Cloudflare challenge: academic.oup.com, jamanetwork.com,
    journals.lww.com, n.neurology.org, and the europepmc.org website (its
    REST API works).
  - OpenGWAS needs a JWT. Open alternatives: GWAS Catalog FTP, FinnGen R10,
    eQTLGen.
  - No email used.
- **Braak 2011:**
  - PDF saved as `data/raw/papers/braak2011.pdf`.
  - The stage-by-age and amyloid-by-age counts are only in SDC 1 and
    SDC 2. Both are open on the publisher's permalinks; saved and extracted.
  - 620 values; 43/43 count checks pass. Two inconsistencies in the source
    are recorded (a printed 34.3% for 34.4%; the Discussion's "52" vs 40
    stage-a cases).
  - Blind second reading from page images: 620/620 agree. `verified` stays
    false.
  - Per-brain amyloid × tau stage × age is not in the paper, so option (b)
    is not possible from it.
- **Module 2:**
  - The real input files are wired in. The pre-fit gate also checks the
    plaque total and the band totals.
  - Prereg data section filled; prereg/config consistency test added.
- **Power study** (100 simulated data sets, first design: whole-band CV):
  - It exposed a structural flaw. A stage seen only in a held-out band gets
    probability ~1e-12 from every model, so the verdict hinges on noise.
  - Fold fits with no start passing the convergence test also stopped the run.
- **Criterion 2 revised** (prereg change log, simulated evidence only):
  - Primary is now CV over brains within age bands.
  - Whole-band CV is kept as sensitivity 7.
  - Unverifiable fold fits are flagged instead of stopping the run.
  - Power study re-run on the new design; results in
    `manuscript/notes/module2_identifiability.md` §6.
- **Module 5:**
  - Validation 4/4.
  - Dry run covers all 5 contrasts (the DHP contrast was missing from the
    simulator).
  - Null calibration of the verdict rule: Knapp–Hartung gives 2.2–3.5% false
    AGAINST per contrast (nominal 2.5%), so it is kept.
- **Module 1:** TwoSampleMR 0.7.12 installed (Posit PM binaries, r-universe).
  - The cross-check is in the pipeline.
  - Estimators and harmonisation agree on simulated data.
  - The GWAS data step is deferred (open, reachable).
- **Modules 3 and 6:** simulated skeletons; Module 6 re-validated 6/6; Module 3 4/5 (published example).
- **Module 2 dry run on the revised design** (simulated, smoke settings), all stages in one pass:
  - primary CV: A −592.0, B −602.7, C −619.9 → supports;
  - whole-band sensitivity: supports;
  - 0 of 30 fold fits unverified.
- **OSF package** `results/osf_package/20261014_dbc4335.zip` (git-ignored).
  - Changed since 20261013_f1fea88: baserate.md, config/baserate.yaml,
    metaanalysis.md, mr.md, config/mr_methods.yaml (listed in
    CHANGES_SINCE_20261013_f1fea88.md inside the package).
  - Not frozen.
