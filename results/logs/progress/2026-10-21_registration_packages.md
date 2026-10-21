# 2026-10-21: verification sheet, registration packages for Modules 1, 3, 4, 6; Module 2 progress

- **Verification sheet** sent: data/extracted/verification/verification_sheet.csv
  (120 Braak counts, then 37 meta rows). `make apply-verification
  TICKS=<sheet> [APPLY=1]` routes the ticks to each file, all or nothing.
- **Registration packages** (plan .md + PDF with text layer + OSF answers;
  each passes src.registration.osf_answers and its consistency test):
  - Module 1 (mr):
    - No target passes its positive control in the ±100 kb window.
    - ADRA2A passes via MetaBrain in the ±1000 kb window only.
    - HMGCR validation passes (Huang 2021).
    - No AD row read.
  - Module 3 (mediation):
    - RADC Standard plan; validation 10/10.
    - braaksc cannot test the LC → cortex step.
    - LC mediation is a conditional branch.
  - Module 4 (dynamics):
    - Validation 12/12; 39/39 machine checks.
    - The pilot found no tipping point (disclosed).
  - Module 6 (network):
    - Validation 6/6.
    - The ADNI DUA forbids AI tools on participant-level data, so the
      researcher runs that step.
- **Module 2** (registered run, resumed from sub-checkpoints after three
  container restarts):
  - Variant A: MLE passes the step rule at h = 0.03125 y; NUTS chain 1 of 4
    done.
  - Expect several more days of wall-clock time.
- **Validation re-run:** metaanalysis re-validated 4/4 after its fingerprint
  changed (MR target files).
