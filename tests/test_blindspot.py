import numpy as np
from cswitch.blindspot import abundance_weighted_rci, is_masked

def test_abundance_weighted_rci():
    # isoform RCI=[-3,+3], 重み=[0.5,0.5] -> 0.0（中間に埋もれる）
    assert abs(abundance_weighted_rci([-3.0, 3.0], [0.5, 0.5]) - 0.0) < 1e-9
    # 重み偏り
    assert abs(abundance_weighted_rci([-3.0, 3.0], [0.9, 0.1]) - (-2.4)) < 1e-9

def test_is_masked_true_when_gene_level_neutral_but_isoforms_split():
    # gene-level が境界付近(|.|<1)だが isoform は両区画にまたがる -> masked
    assert is_masked(gene_rci=0.1, isoform_rcis=[-2.5, 2.0], neutral=1.0) is True

def test_is_masked_false_when_concordant():
    assert is_masked(gene_rci=-2.0, isoform_rcis=[-2.5, -1.8], neutral=1.0) is False
