import numpy as np

def abundance_weighted_rci(rcis, weights):
    rcis = np.asarray(rcis, float); w = np.asarray(weights, float)
    return float(np.sum(rcis * w) / np.sum(w))

def is_masked(gene_rci, isoform_rcis, neutral=1.0):
    """gene-level が中立帯(|gene_rci|<neutral)に居るのに、isoform が両区画にまたがる。"""
    arr = np.asarray(isoform_rcis, float)
    spans_both = (arr.max() > 0) and (arr.min() < 0)
    return bool(abs(gene_rci) < neutral and spans_both)
