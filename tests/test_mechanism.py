import numpy as np, pandas as pd
from cswitch.mechanism import alu_content_diff, length_matched_compare

def test_alu_content_diff():
    # isoform A: Alu塩基=200/長1000=0.2, B: 0/長500=0.0 -> diff=0.2
    assert abs(alu_content_diff(0.2, 0.0) - 0.2) < 1e-9

def test_length_matched_compare_returns_pvalue():
    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "group": ["cs"]*30 + ["concordant"]*30,
        "alu_diff": np.r_[rng.normal(0.15,0.05,30), rng.normal(0.05,0.05,30)],
        "len_diff": np.r_[rng.normal(500,100,30), rng.normal(500,100,30)],
    })
    res = length_matched_compare(df, value="alu_diff", group="group", covariate="len_diff")
    assert "p_value" in res and 0.0 <= res["p_value"] <= 1.0
    assert res["median_cs"] > res["median_concordant"]
