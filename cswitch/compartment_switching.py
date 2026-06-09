import numpy as np, pandas as pd


def delta_rci(rci_a, rci_b):
    return float(rci_b - rci_a)


def classify_compartment_switch(rci_a, rci_b, threshold=1.0):
    """符号反転（反対区画）かつ |ΔRCI| >= threshold。"""
    sign_flip = (np.sign(rci_a) != np.sign(rci_b)) and (rci_a != 0) and (rci_b != 0)
    return bool(sign_flip and abs(rci_b - rci_a) >= threshold)


def permutation_null_fraction(observed_fraction, iso_rci, n_pairs, threshold, B=10000, seed=0):
    """同一遺伝子内ランダム isoform ペアで compartment-switch 割合の帰無分布を作る。"""
    rng = np.random.default_rng(seed)
    by_gene = {g: d["rci"].values for g, d in iso_rci.groupby("gene_id") if len(d) >= 2}
    genes = list(by_gene.keys())
    null = np.empty(B)
    for b in range(B):
        hits = 0
        for _ in range(n_pairs):
            g = genes[rng.integers(len(genes))]
            vals = by_gene[g]
            i, j = rng.choice(len(vals), size=2, replace=False)
            if classify_compartment_switch(vals[i], vals[j], threshold):
                hits += 1
        null[b] = hits / n_pairs
    p = float((np.sum(null >= observed_fraction) + 1) / (B + 1))
    return observed_fraction, p, null


def permutation_null_gene_fraction(observed_fraction, gene_events, gene_pool,
                                   threshold, B=10000, seed=0):
    """Matched null for the gene-level compartment-switcher fraction.

    Reproduces the OBSERVED statistic exactly: the fraction of genes with at least
    one compartment-switch event. For each gene, its observed dominant-isoform pairs
    are replaced by random within-gene isoform pairs (drawn from the gene's localized
    isoforms), preserving the per-gene number of events; a gene is a hit if ANY of its
    drawn pairs qualifies. One-sided p = P(null fraction >= observed).

    gene_events : dict gene_id -> number of switch events for that gene
    gene_pool   : dict gene_id -> 1D array of localized isoform RCIs (>= 2 values)
    """
    rng = np.random.default_rng(seed)
    genes = [g for g in gene_events
             if g in gene_pool and len(gene_pool[g]) >= 2]
    null = np.empty(B)
    for b in range(B):
        hits = 0
        for g in genes:
            pool = gene_pool[g]
            n = len(pool)
            gene_hit = False
            for _ in range(gene_events[g]):
                i, j = rng.choice(n, size=2, replace=False)
                if classify_compartment_switch(pool[i], pool[j], threshold):
                    gene_hit = True
                    break
            if gene_hit:
                hits += 1
        null[b] = hits / len(genes)
    p = float((np.sum(null >= observed_fraction) + 1) / (B + 1))
    return observed_fraction, p, null
