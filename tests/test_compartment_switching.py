import numpy as np, pandas as pd
from cswitch.compartment_switching import (
    delta_rci, classify_compartment_switch,
    permutation_null_fraction, permutation_null_gene_fraction,
)

def test_delta_rci():
    assert abs(delta_rci(rci_a=-2.0, rci_b=1.5) - 3.5) < 1e-9

def test_classify_sign_flip_and_threshold():
    # 符号反転 & |Δ|>=1 -> compartment switch
    assert classify_compartment_switch(-2.0, 1.5, threshold=1.0) is True
    # 同符号（両方核）-> False
    assert classify_compartment_switch(-2.0, -0.5, threshold=1.0) is False
    # 符号反転だが |Δ|<threshold -> False
    assert classify_compartment_switch(-0.2, 0.2, threshold=1.0) is False

def test_permutation_null_fraction_reproducible():
    iso = pd.DataFrame({
        "gene_id": ["G1"]*4 + ["G2"]*3,
        "transcript_id": ["t1","t2","t3","t4","u1","u2","u3"],
        "rci": [-3.0, 2.0, -1.0, 1.0, -2.0, -1.5, 0.5],
    })
    obs, p, null = permutation_null_fraction(
        observed_fraction=0.5, iso_rci=iso, n_pairs=2, threshold=1.0, B=200, seed=42)
    assert 0.0 <= p <= 1.0
    assert len(null) == 200
    _, p2, _ = permutation_null_fraction(0.5, iso, 2, 1.0, B=200, seed=42)
    assert p == p2

def test_permutation_null_gene_fraction_matched():
    # G1 pool can cross (pairs (-3,3),(-3,1) cross; (3,1) does not -> 2/3 cross)
    # G2 pool [-2,-1] same compartment -> never crosses
    gene_pool = {"G1": np.array([-3.0, 3.0, 1.0]), "G2": np.array([-2.0, -1.0])}
    gene_events = {"G1": 1, "G2": 1}
    obs, p, null = permutation_null_gene_fraction(
        0.5, gene_events, gene_pool, threshold=1.0, B=500, seed=1)
    assert len(null) == 500
    assert 0.0 <= p <= 1.0
    # only G1 can ever be a hit -> fraction is 0 or 0.5, never above 0.5
    assert null.max() <= 0.5 + 1e-9
    # G2 never hits; mean fraction ~ (2/3)/2
    assert 0.2 < null.mean() < 0.45
    # deterministic with same seed
    _, p2, _ = permutation_null_gene_fraction(0.5, gene_events, gene_pool, 1.0, B=500, seed=1)
    assert p == p2
