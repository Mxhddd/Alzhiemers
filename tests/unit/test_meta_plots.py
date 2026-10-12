import numpy as np
import pandas as pd

from src.metaanalysis.plots import forest_figure, triangulation_figure


def test_triangulation_figure_keeps_missing_estimates_as_gaps(tmp_path):
    tri = pd.DataFrame([
        {"target": "ADRB1", "prediction_lower_activity": "protective", "obs_est": -0.2, "obs_ci_lower": -0.3,
         "obs_ci_upper": -0.1, "mr_est": -0.1, "mr_ci_lower": -0.15, "mr_ci_upper": -0.05,
         "obs_vs_prediction": "agree", "mr_vs_prediction": "agree", "obs_vs_mr": "same sign"},
        {"target": "LGMN", "prediction_lower_activity": "protective", "obs_est": np.nan, "obs_ci_lower": np.nan,
         "obs_ci_upper": np.nan, "mr_est": np.nan, "mr_ci_lower": np.nan, "mr_ci_upper": np.nan,
         "obs_vs_prediction": "missing", "mr_vs_prediction": "missing", "obs_vs_mr": "missing"},
    ])
    fig = triangulation_figure(tri)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert "no observational estimate" in texts and "no MR estimate" in texts
    assert [t.get_text() for t in fig.axes[0].get_yticklabels()] == ["LGMN", "ADRB1"]


def test_forest_handles_unpooled_rows():
    d = pd.DataFrame([{"group": "a", "k": 3, "est": -0.1, "ci_lower": -0.2, "ci_upper": 0.0, "pi_lower": -0.4,
                       "pi_upper": 0.2, "status": "pooled"},
                      {"group": "b", "k": 0, "est": np.nan, "ci_lower": np.nan, "ci_upper": np.nan,
                       "pi_lower": np.nan, "pi_upper": np.nan, "status": "not pooled"}])
    assert len(forest_figure(d).axes) == 1
