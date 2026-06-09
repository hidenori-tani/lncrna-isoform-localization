import numpy as np, pandas as pd

def alu_content_diff(frac_a, frac_b):
    return float(abs(frac_a - frac_b))

def length_matched_compare(df, value, group, covariate, n_bins=5, B=10000, seed=0):
    """Length-matched group comparison via WITHIN-BIN label permutation.

    Bin by quantiles of `covariate` (e.g. transcript-length difference); the statistic
    is median(cs) - median(concordant) of `value`. The null permutes the cs/concordant
    labels WITHIN each length bin (so the length distribution is matched), giving a
    one-sided p for cs > concordant. This genuinely uses the bins (unlike a plain,
    unstratified Mann-Whitney).
    """
    d = df.dropna(subset=[value, covariate]).copy()
    nb = min(n_bins, d[covariate].nunique())
    d["bin"] = pd.qcut(d[covariate], q=nb, duplicates="drop") if nb >= 2 else 0
    vals = d[value].to_numpy(float)
    is_cs = (d[group].to_numpy() == "cs")
    bins = d["bin"].astype(str).to_numpy()

    def stat(mask):
        cs_v = vals[mask]; cc_v = vals[~mask]
        if len(cs_v) == 0 or len(cc_v) == 0:
            return 0.0
        return float(np.median(cs_v) - np.median(cc_v))

    obs = stat(is_cs)
    rng = np.random.default_rng(seed)
    idx_by_bin = [np.where(bins == b)[0] for b in pd.unique(bins)]
    null = np.empty(B)
    for k in range(B):
        perm = np.empty(len(vals), dtype=bool)
        for idx in idx_by_bin:
            perm[idx] = rng.permutation(is_cs[idx])
        null[k] = stat(perm)
    p = float((np.sum(null >= obs) + 1) / (B + 1))
    return {"p_value": p, "median_cs": float(np.median(vals[is_cs])),
            "median_concordant": float(np.median(vals[~is_cs])), "obs_stat": obs}
