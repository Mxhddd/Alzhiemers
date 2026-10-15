# 2026-10-15 — OSF check, freezes, first registered runs

- **OSF check.**
  - xb7fy = baserate and 9kb4w = metaanalysis, matched by title. Both are
    public and not withdrawn, and each has a document attached.
  - The registered files are text-less PDFs, so a byte-for-byte match with
    the Markdown is impossible. The PDFs are kept byte-identical in
    prereg/registered/ (sha256 = OSF's own hash). prereg/*.md are now
    transcriptions of the registered text.
  - All registered numbers and rules match config. Tests: 195 Python, 218 R.
- **Frozen:** baserate and metaanalysis only.
- **Module 5** (real run; UNVERIFIED inputs):
  - The positive control was recovered: 1.62 (1.26–2.08), k = 7.
  - β-blockers: inconclusive (two studies disagree, I² 91%).
  - α1-blockers: inconclusive (not robust; supports under S2 and S4).
  - DHP: supports, 0.82 (0.72–0.93), one study, robust.
  - Atomoxetine: no data.
  - No AGAINST. Overall: neither against nor in favour by the registered
    rule. No departures.
- **Module 2:** running.
  - The registered step rule needs h = 0.03125 y on the real data (bound
    2.3e-2 at 0.25 y; the Radau check passes at every step). That is about
    two days of compute.
  - The first attempt was killed by the 2 h tool limit before a checkpoint
    was written. It was restarted unchanged as a detached process.
- **Module 1 data step:** done; no AD outcome rows were read.
  - validate-mr fails: the targets are incomplete and need the researcher's
    values.
  - All 4 positive controls that could run fail (ADRB1, ADRA2A, DBH: not
    significant; CACNA1D: wrong sign). 4 genes are absent from eQTLGen.
    See manuscript/notes/module1_data_step.md.
- **Verification checklist + tick script:** data/extracted/verification/;
  make verification-checklist / apply-verification.
- **Departures:** manuscript/notes/deviations.md.
