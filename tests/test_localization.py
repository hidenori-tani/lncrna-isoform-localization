import numpy as np, pandas as pd
from cswitch.localization import cn_rci, consensus_rci

def test_cn_rci_basic():
    # 細胞質 TPM=8, 核 TPM=2 -> log2(8/2)=2.0
    assert abs(cn_rci(8.0, 2.0) - 2.0) < 1e-9
    # 核優位: 細胞質2,核8 -> -2.0
    assert abs(cn_rci(2.0, 8.0) + 2.0) < 1e-9

def test_cn_rci_pseudocount_guards_zero():
    v = cn_rci(0.0, 10.0, pseudocount=0.01)
    assert np.isfinite(v) and v < 0

def test_cn_rci_below_min_tpm_returns_nan():
    assert np.isnan(cn_rci(0.001, 0.002, min_tpm=0.1))

def test_consensus_rci_median_and_agreement():
    df = pd.DataFrame({
        "transcript_id": ["T1","T1","T1","T2","T2"],
        "cell_line": ["A","B","C","A","B"],
        "rci": [2.0, 2.4, 1.6, -3.0, -2.0],
    })
    out = consensus_rci(df).set_index("transcript_id")
    assert abs(out.loc["T1","rci_consensus"] - 2.0) < 1e-9   # median
    assert out.loc["T1","n_cell_lines"] == 3
    assert abs(out.loc["T1","sign_agreement"] - 1.0) < 1e-9
