"""Task 10: compartment-switching main analysis.

- ΔRCI per switch event; compartment-switch classification (sign flip + |ΔRCI|>=thr)
- Headline: fraction of localized switch genes that are compartment-switchers
- Permutation null: random within-gene isoform pairs (from ALL localized isoforms of
  the analysis genes, not just the dominant ones) -> tests whether the observed
  dominant-pair switching crosses compartments more than chance.
- Threshold sensitivity (0.5/1.0/1.5/2.0).
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd

from cswitch.compartment_switching import (
    delta_rci,
    classify_compartment_switch,
    permutation_null_gene_fraction,
)
from cswitch.keynumbers import set_number

SWL = os.path.join(ROOT, "data/processed/switch_with_localization.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
OUT = os.path.join(ROOT, "data/processed/compartment_calls.tsv")

PRIMARY_THR = 1.0
THRESHOLDS = [0.5, 1.0, 1.5, 2.0]


def gene_fraction(df, thr):
    cs = [classify_compartment_switch(a, b, thr) for a, b in zip(df.rci_a, df.rci_b)]
    d = df.assign(cs=cs)
    return d, d.groupby("gene_id")["cs"].any()


def main():
    df = pd.read_csv(SWL, sep="\t")
    df["delta_rci"] = [delta_rci(a, b) for a, b in zip(df.rci_a, df.rci_b)]

    # primary threshold
    d, gene_cs = gene_fraction(df, PRIMARY_THR)
    d.to_csv(OUT, sep="\t", index=False)
    frac = float(gene_cs.mean())
    n_cs_genes = int(gene_cs.sum())
    n_genes = int(gene_cs.size)
    n_cs_events = int(d.cs.sum())
    n_events = int(len(d))

    # Matched null (Codex review fix): reproduce the OBSERVED statistic exactly =
    # fraction of genes with >=1 compartment-switch event. For each gene, replace its
    # dominant-isoform pairs with random within-gene isoform pairs, preserving the
    # per-gene number of events; a gene is a hit if ANY drawn pair qualifies.
    # gene keys are version-stripped to match isoform_localization.
    d["gid"] = d.gene_id.map(lambda x: str(x).split(".")[0])
    gene_events = d.groupby("gid").size().to_dict()
    loc = pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"])
    loc = loc[loc.gene_id.isin(set(gene_events))]
    gene_pool = {g: grp.rci_consensus.to_numpy()
                 for g, grp in loc.groupby("gene_id") if len(grp) >= 2}

    _, p, null = permutation_null_gene_fraction(
        observed_fraction=frac, gene_events=gene_events, gene_pool=gene_pool,
        threshold=PRIMARY_THR, B=10000, seed=42,
    )

    set_number("compartment_switcher_gene_fraction", frac)
    set_number("compartment_switcher_n_genes", n_cs_genes)
    set_number("compartment_switcher_n_genes_total", n_genes)
    set_number("compartment_switch_events", n_cs_events)
    set_number("compartment_switch_events_total", n_events)
    set_number("compartment_switcher_permutation_p", p)
    set_number("compartment_null_pool_genes", int(len(gene_pool)))
    set_number("compartment_null_mean", float(null.mean()))
    pd.DataFrame({"null_fraction": null}).to_csv(
        os.path.join(ROOT, "data/processed/null_distribution.tsv"), sep="\t", index=False)

    print(f"threshold={PRIMARY_THR}")
    print(f"  compartment-switch genes = {n_cs_genes}/{n_genes} (frac={frac:.3f})")
    print(f"  compartment-switch events = {n_cs_events}/{n_events}")
    print(f"  permutation p = {p:.4g} (null pool genes={len(gene_pool)}, "
          f"null mean={null.mean():.3f})")

    # sensitivity
    sens = {}
    for thr in THRESHOLDS:
        _, gc = gene_fraction(df, thr)
        sens[str(thr)] = float(gc.mean())
    set_number("compartment_switcher_fraction_by_threshold", sens)
    print("sensitivity (gene fraction):", sens)


if __name__ == "__main__":
    main()
