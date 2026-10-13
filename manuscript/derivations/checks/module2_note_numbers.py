"""Machine check of the numbers quoted in manuscript/notes/module2_identifiability.md
that are arithmetic rather than simulation output.

    uv run python manuscript/derivations/checks/module2_note_numbers.py
"""
import sys

from scipy.stats import norm

results = []


def check(name, ok, detail=""):
    results.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{(': ' + detail) if detail else ''}")


se0 = 0.388                                   # SE of log HR at n = 2,300, no age slope (simulation output)
thresh = 0.322 - 1.959963984540054 * se0      # estimate needed for the upper 95% limit to fall below 0.322
p_against = norm.cdf(thresh / se0)
check("'against' needs the estimate below -0.44", abs(thresh - (-0.44)) < 0.005, f"{thresh:.4f}")
check("which happens ~13% of the time under the Wald approximation", abs(p_against - 0.13) < 0.005, f"{p_against:.4f}")
check("dry-run in-sample gap A vs B = 9.1", abs((-166.55) - (-175.63) - 9.1) < 0.05, f"{-166.55 + 175.63:.2f}")
print(f"\n{sum(results)}/{len(results)} checks passed")
sys.exit(0 if all(results) else 1)
