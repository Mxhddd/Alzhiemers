# Module 5 — hand-computable three-study example

> **Status: derived by Claude; machine-checked (exact rational arithmetic in
> sympy, `checks/module5_hand_three_study.py`), NOT human-checked.** The
> brief planned for the researcher to supply this arithmetic. Under the
> 2026-10-13 instruction it was done by Claude. Result that depends on it:
> the validation target `tests/validation/targets/meta_hand_three_study.yaml`.

Inputs are test settings chosen to be exact rationals, not published values:

| study | y (log ratio) | v (variance) | SE | w = 1/v |
|---|---|---|---|---|
| 1 | 1/5 = 0.2 | 1/25 | 0.2 | 25 |
| 2 | 1/2 = 0.5 | 1/16 | 0.25 | 16 |
| 3 | −1/10 = −0.1 | 1/100 | 0.1 | 100 |

**Fixed effect (inverse variance).**
- Σw = 141.
- μ̂_FE = (25·0.2 + 16·0.5 + 100·(−0.1)) / 141 = (5 + 8 − 10)/141 = 3/141 = **1/47 ≈ 0.0212765957**.
- SE_FE = 1/√141 ≈ **0.0842151921**.

**Cochran's Q.** Q = Σ w (y − μ̂)², with y − μ̂ equal to 42/235, 45/94 and −57/470:
- 25·(42/235)² = 25·1764/55225 = 1764/2209;
- 16·(45/94)² = 16·2025/8836 = 8100/2209;
- 100·(57/470)² = 100·3249/220900 = 3249/2209.
- Sum: Q = 13113/2209 = **279/47 ≈ 5.9361702128**.

**DerSimonian–Laird τ².**
- C = Σw − Σw²/Σw = 141 − 10881/141 = 9000/141.
- τ² = (Q − (k − 1))/C = (279/47 − 2)·141/9000 = (185/47)·(141/9000) = 555/9000 = **37/600 ≈ 0.0616666667**.

**Random effects (DL).**
- Weights w* = 1/(v + τ²): v + τ² equals 61/600, 149/1200 and 43/600, so w* = 600/61, 1200/149 and 600/43.
- μ̂_RE = Σw*y/Σw* = **1997/13828 ≈ 0.1444171247**.
- SE_RE = (Σw*)^−1/2 = √1351088939/207420 ≈ **0.1772112706**.

The pipeline (metafor `rma`, methods "FE" and "DL") must reproduce all six
numbers to 1e-9. That is far below any rounding: the expected values are
exact.

Machine check (`checks/module5_hand_three_study.py`) output:

    fe_estimate     exact = 1/47   ~ 0.0212765957446809
    fe_se           exact = sqrt(141)/141   ~ 0.0842151921066519
    q               exact = 279/47   ~ 5.93617021276596
    tau2_dl         exact = 37/600   ~ 0.0616666666666667
    re_dl_estimate  exact = 1997/13828   ~ 0.144417124674573
    re_dl_se        exact = sqrt(1351088939)/207420   ~ 0.177211270550306
