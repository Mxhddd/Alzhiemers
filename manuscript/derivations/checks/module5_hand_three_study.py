"""Exact arithmetic for the Module 5 hand-computable three-study example
(manuscript/derivations/module5_hand_three_study.md). Inputs are TEST
SETTINGS chosen to be exact rationals, not published values. Every quantity
is computed in exact rational arithmetic (sympy), then printed to 15
significant digits; the square roots are the only irrational steps.

    uv run python manuscript/derivations/checks/module5_hand_three_study.py
"""
from sympy import Integer, Rational, sqrt

y = [Rational(1, 5), Rational(1, 2), Rational(-1, 10)]      # log ratios
v = [Rational(1, 25), Rational(1, 16), Rational(1, 100)]    # sampling variances (SE 0.2, 0.25, 0.1)
k = len(y)

w = [1 / vi for vi in v]                                    # 25, 16, 100
W = sum(w)                                                  # 141
fe = sum(wi * yi for wi, yi in zip(w, y)) / W               # 3/141 = 1/47
fe_se = 1 / sqrt(W)
Q = sum(wi * (yi - fe) ** 2 for wi, yi in zip(w, y))
C = W - sum(wi ** 2 for wi in w) / W
tau2 = max(Integer(0), (Q - (k - 1)) / C)
ws = [1 / (vi + tau2) for vi in v]
re = sum(wi * yi for wi, yi in zip(ws, y)) / sum(ws)
re_se = 1 / sqrt(sum(ws))

out = {"fe_estimate": fe, "fe_se": fe_se, "q": Q, "tau2_dl": tau2, "re_dl_estimate": re, "re_dl_se": re_se}
for name, val in out.items():
    print(f"{name:15s} exact = {val}   ~ {val.evalf(15)}")
