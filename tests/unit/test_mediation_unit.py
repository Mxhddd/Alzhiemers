import numpy as np
import pytest

from src.mediation.mediation import estimate, killer2_verdict


def test_bad_inputs_raise():
    T = np.array([0.0, 1, 0, 1]); M = np.array([1.0, 2, 1, 2]); Y = np.array([1.0, 2, 3, 4])
    with pytest.raises(ValueError, match="non-finite"):
        estimate(T, M, np.array([1.0, np.nan, 3, 4]), None)
    with pytest.raises(ValueError, match="no variation"):
        estimate(np.zeros(4), M, Y, None)
    with pytest.raises(ValueError, match="rho"):
        estimate(T, M + np.array([0.1, -0.1, 0.2, 0.0]), Y, None).acme_at(np.array([1.0]))


@pytest.mark.parametrize("acme,lo,hi,expect", [
    (0.2, 0.1, 0.3, "NOT met"),
    (-0.2, -0.3, -0.1, "OPPOSITE"),
    (0.05, -0.02, 0.12, "includes 0"),
])
def test_killer2_verdicts(acme, lo, hi, expect):
    s = {"estimates": {"acme": {"estimate": acme, "ci_lower": lo, "ci_upper": hi}}}
    assert expect in killer2_verdict(s, predicted_sign=1)


def test_real_path_fails_loudly_and_names_the_provider():
    from src.common.io import MissingInputError
    from src.common.config import load_config
    from src.common.paths import CONFIG_DIR
    from src.mediation.run_mediation import load_cohort

    class Run:
        def drop(self, **kw): raise AssertionError("no drops expected")
    with pytest.raises(MissingInputError, match="no access yet"):
        load_cohort(load_config(CONFIG_DIR / "mediation.yaml"), Run())
