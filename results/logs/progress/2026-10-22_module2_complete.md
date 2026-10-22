# 2026-10-22: Module 2 registered run complete (provisional: unverified inputs)

- **Criterion 2** (5-fold CV over brains within age bands; margin 2):
  - **AGAINST**.
  - Held-out log-likelihood: A −546.72, B −543.63, C −555.73.
  - B beats A by 3.09. 0/15 unverified folds.
- **Falsification** (§6, uncertainty-aware): PASS for A, B and C.
  - The point-estimate rule FAILS for all three (85+ band).
  - All variants under-predict incidence from age 75 on, but stay within
    the factor-2 tolerance.
- **Overall (§2 rule):** counts AGAINST the LC-first model, because
  criterion 2 is "against". Not called robust: criterion 2 was not
  re-evaluated under analyses 1–5 (deviations.md note 6).
- **Sensitivity analyses 1–6:** falsification PASS in every one.
- **Sensitivity 7** (whole-band CV): inconclusive.
  - A −218.20, B −218.87, C −237.63.
  - B and C each had 1/5 unverified fold.
- **Descriptive amyloid HR** (variant A, NUTS):
  - log HR −0.03 (95% HDI −2.03 to 1.78).
  - Labelled "not identified: estimate on the parameter bound": the MLE
    sits at the −5 bound, and the posterior ≈ the N(0,1) prior.
  - All sensitivity estimates are on the bound.
- **Diagnostics:** all three posteriors clean (R-hat 1.00, 0 divergences).
  - Steps: A 0.03125 y, B 0.0625 y, C 0.125 y.
- **Run history:** five container restarts; resumed from checkpoints each
  time.
