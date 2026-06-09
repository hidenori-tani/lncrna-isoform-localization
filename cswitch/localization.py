import numpy as np
import pandas as pd


def cn_rci(cyto_tpm, nuc_tpm, pseudocount=0.01, min_tpm=0.1):
    """CN-RCI = log2(cytoplasm/nucleus). 両分画とも min_tpm 未満なら NaN。
    pseudocount はゼロ除算防止のためにのみ加算する（値が 0 のとき）。"""
    if (cyto_tpm < min_tpm) and (nuc_tpm < min_tpm):
        return np.nan
    c = cyto_tpm if cyto_tpm != 0.0 else pseudocount
    n = nuc_tpm if nuc_tpm != 0.0 else pseudocount
    return float(np.log2(c / n))


def consensus_rci(df):
    """df: transcript_id, cell_line, rci -> 細胞株間 median + 一致度。"""
    def agg(g):
        vals = g["rci"].dropna().values
        if len(vals) == 0:
            return pd.Series({
                "rci_consensus": np.nan,
                "n_cell_lines": 0,
                "sign_agreement": np.nan,
            })
        signs = np.sign(vals)
        signs = signs[signs != 0]
        if len(signs):
            agree = float(max(np.mean(signs > 0), np.mean(signs < 0)))
        else:
            agree = np.nan
        return pd.Series({
            "rci_consensus": float(np.median(vals)),
            "n_cell_lines": int(len(vals)),
            "sign_agreement": agree,
        })

    return (
        df.groupby("transcript_id")[["cell_line", "rci"]]
        .apply(agg)
        .reset_index()
    )
