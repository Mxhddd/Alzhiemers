import math

import pandas as pd
import pytest

from src.mr.forest_plot import draw, forest_rows


def _summary():
    return pd.DataFrame([
        {"gene": "ADRB1", "exposure_source": "eqtlgen", "predicted_lower_activity": "protective",
         "method": "ivw", "nsnp": 3, "b": 0.1, "ci_lower": 0.02, "ci_upper": 0.18, "status": "estimated",
         "is_primary": True, "analysis": "main", "concordance": "concordant"},
        {"gene": "ADRB1", "exposure_source": "eqtlgen", "predicted_lower_activity": "protective",
         "method": "mr_egger", "nsnp": 3, "b": 0.3, "ci_lower": -0.5, "ci_upper": 1.1, "status": "estimated",
         "is_primary": False, "analysis": "main", "concordance": None},
        {"gene": "SLC6A2", "exposure_source": "eqtlgen", "predicted_lower_activity": "protective",
         "method": None, "nsnp": None, "b": None, "ci_lower": None, "ci_upper": None,
         "status": "excluded: no positive control specified", "is_primary": None, "analysis": None,
         "concordance": None},
    ])


def test_one_row_per_target_and_excluded_targets_kept_as_gaps():
    rows = forest_rows(_summary())
    assert list(rows["label"]) == ["ADRB1 [eqtlgen]", "SLC6A2 [eqtlgen]"]
    assert rows.loc[0, "b"] == 0.1
    assert math.isnan(rows.loc[1, "b"]) and "no positive control" in rows.loc[1, "note"]


def test_ambiguous_summary_is_refused():
    bad = _summary()
    bad.loc[1, "is_primary"] = True   # two primary estimates for one target
    with pytest.raises(ValueError, match="expected one primary"):
        forest_rows(bad)


def test_draw_handles_gaps():
    fig = draw(forest_rows(_summary()))
    assert len(fig.axes) == 1
