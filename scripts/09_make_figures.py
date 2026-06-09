"""Task 16: figures (honest thesis).

Fig1  isoform-resolved localization landscape of switch-gene isoforms
Fig2  compartment crossing is common but NOT switching-enriched (vs permutation null) + threshold sensitivity
Fig3  gene-level localization blind spot (core): gene-level masks isoform compartment split + robustness curve
Fig4  Alu does NOT explain compartment differences (honest negative)
Fig5  robust case studies

All numbers from data/processed/* and key_numbers.json (no hardcoded results).
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = os.path.join(ROOT, "figures")
PROC = os.path.join(ROOT, "data/processed")
KN = json.load(open(os.path.join(ROOT, "key_numbers.json")))

os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 300, "font.size": 10, "axes.spines.top": False, "axes.spines.right": False})

loc = pd.read_csv(os.path.join(PROC, "isoform_localization.tsv"), sep="\t").dropna(subset=["rci_consensus"])
calls = pd.read_csv(os.path.join(PROC, "compartment_calls.tsv"), sep="\t")
bs = pd.read_csv(os.path.join(PROC, "blindspot.tsv"), sep="\t")
pairs = pd.read_csv(os.path.join(PROC, "mechanism_pairs.tsv"), sep="\t")
cases = pd.read_csv(os.path.join(PROC, "case_studies.tsv"), sep="\t")

NUC = "#2c6fbb"  # nuclear (negative RCI)
CYT = "#d1495b"  # cytoplasmic (positive RCI)


def fig1():
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ana = set(calls.gene_id.map(lambda x: str(x).split(".")[0]))
    r = loc[loc.gene_id.isin(ana)].rci_consensus
    ax[0].hist(r, bins=40, color="#7a7a7a")
    ax[0].axvline(0, color="k", lw=1, ls="--")
    ax[0].set_xlabel("isoform CN-RCI  (log2 cyto/nuc)")
    ax[0].set_ylabel("isoforms")
    ax[0].set_title(f"A  Isoform localization\n(n={len(r)} isoforms, {loc[loc.gene_id.isin(ana)].gene_id.nunique()} switch genes)")
    ax[0].text(0.02, 0.95, "nuclear", color=NUC, transform=ax[0].transAxes, va="top")
    ax[0].text(0.98, 0.95, "cytoplasmic", color=CYT, transform=ax[0].transAxes, va="top", ha="right")
    # per-gene isoform RCI range
    g = loc[loc.gene_id.isin(ana)].groupby("gene_id").rci_consensus.agg(["min", "max", "count"])
    g = g[g["count"] >= 2].sort_values("min")
    y = np.arange(len(g))
    ax[1].hlines(y, g["min"], g["max"], color="#bbbbbb", lw=0.6)
    ax[1].scatter(g["min"], y, s=6, color=NUC)
    ax[1].scatter(g["max"], y, s=6, color=CYT)
    ax[1].axvline(0, color="k", lw=1, ls="--")
    ax[1].set_xlabel("isoform CN-RCI range per gene")
    ax[1].set_ylabel("switch genes (sorted)")
    ax[1].set_title(f"B  Within-gene isoform spread\n({len(g)} genes, >=2 localized isoforms)")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig1.png")); plt.close(fig)


def fig2():
    # matched null distribution saved by scripts/05_compartment_analysis.py
    null = pd.read_csv(os.path.join(PROC, "null_distribution.tsv"), sep="\t")["null_fraction"].to_numpy()
    obs = KN["compartment_switcher_gene_fraction"]
    p = KN["compartment_switcher_permutation_p"]

    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ax[0].hist(null, bins=40, color="#bbbbbb", label="within-gene null (all isoforms)")
    ax[0].axvline(obs, color=CYT, lw=2, label=f"observed = {obs:.2f}")
    ax[0].set_xlabel("compartment-switcher gene fraction")
    ax[0].set_ylabel("permutations")
    ax[0].set_title(f"A  Compartment crossing vs within-gene null\n(observed {obs:.2f} < null mean {null.mean():.2f}; one-sided p = {p:.2f})")
    ax[0].legend(fontsize=8)
    sens = KN["compartment_switcher_fraction_by_threshold"]
    ks = sorted(sens, key=float)
    ax[1].plot([float(k) for k in ks], [sens[k] for k in ks], "o-", color="#444")
    ax[1].set_ylim(0, max(sens.values()) * 1.3)
    ax[1].set_xlabel("|ΔRCI| threshold")
    ax[1].set_ylabel("compartment-switcher gene fraction")
    ax[1].set_title("B  Threshold sensitivity")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2.png")); plt.close(fig)


def fig3():
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    masked = bs[bs.masked]; other = bs[~bs.masked]
    ax[0].hlines(other.index, other.iso_rci_min, other.iso_rci_max, color="#dddddd", lw=0.6)
    ax[0].hlines(masked.index, masked.iso_rci_min, masked.iso_rci_max, color="#f0c46b", lw=0.8)
    ax[0].scatter(bs.gene_rci, bs.index, s=10, color="k", label="gene-level RCI", zorder=3)
    ax[0].axvline(0, color="k", lw=1, ls="--")
    ax[0].axvspan(-1, 1, color="#cccccc", alpha=0.3)
    ax[0].set_xlabel("CN-RCI")
    ax[0].set_ylabel("switch genes")
    ax[0].set_title("A  Gene-level RCI (dots) vs isoform range (bars)\nyellow = masked (gene neutral, isoforms split)")
    ax[0].legend(fontsize=8, loc="lower right")
    # robustness curve
    labels = ["all", "agree>=0.8\nncl>=3", "agree>=1.0\nncl>=5"]
    vals = [KN["blindspot_genes_span_both_compartments"] / KN["blindspot_genes_tested"],
            KN["blindspot_span_both_fraction_robust_a08_n3"],
            KN["blindspot_span_both_fraction_robust_a10_n5"]]
    ax[1].bar(labels, [v * 100 for v in vals], color=["#9aa0a6", "#6b8fb5", "#2c6fbb"])
    ax[1].set_ylabel("% genes with isoforms in BOTH compartments")
    ax[1].set_title(f"B  Robust to localization confidence\n(masked overall = {KN['blindspot_masked_fraction']*100:.0f}% of genes)")
    for i, v in enumerate(vals):
        ax[1].text(i, v * 100 + 1, f"{v*100:.0f}%", ha="center", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3.png")); plt.close(fig)


def fig4():
    fig, ax = plt.subplots(figsize=(4.2, 3.8))
    data = [pairs[pairs.group == "concordant"].alu_diff.dropna(),
            pairs[pairs.group == "cs"].alu_diff.dropna()]
    ax.boxplot(data, tick_labels=["concordant", "compartment\nswitch"], showfliers=False)
    ax.set_ylabel("|Alu content difference| between isoforms")
    ax.set_title(f"Alu difference: cs vs concordant (all pairs)\n(p = {KN['mechanism_alu_diff_p']:.2f} overall; p = {KN['mechanism_alu_diff_p_nonzero']:.2f} among Alu-bearing, n={KN['mechanism_alu_nonzero_n']})")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4.png")); plt.close(fig)


def fig5():
    sel = cases.head(3)
    fig, axes = plt.subplots(1, len(sel), figsize=(3.2 * len(sel), 3.6), squeeze=False)
    for ax, (_, row) in zip(axes[0], sel.iterrows()):
        gid = row.gid
        giso = loc[loc.gene_id == gid]
        ax.scatter(giso.rci_consensus, np.zeros(len(giso)) + 0.5, s=30, color="#999", zorder=2)
        a = str(row.dom_isoform_A).split(".")[0]; b = str(row.dom_isoform_B).split(".")[0]
        ax.scatter([row.rci_a], [0.5], s=90, color=(NUC if row.rci_a < 0 else CYT), zorder=3, edgecolor="k")
        ax.scatter([row.rci_b], [0.5], s=90, color=(NUC if row.rci_b < 0 else CYT), zorder=3, edgecolor="k")
        ax.annotate(row.tissueA.split(" - ")[-1][:14], (row.rci_a, 0.5), (row.rci_a, 0.72),
                    ha="center", fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6))
        ax.annotate(row.tissueB.split(" - ")[-1][:14], (row.rci_b, 0.5), (row.rci_b, 0.28),
                    ha="center", fontsize=7, arrowprops=dict(arrowstyle="-", lw=0.6))
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_xlim(min(giso.rci_consensus.min(), -8) - 0.5, max(giso.rci_consensus.max(), 8) + 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel("isoform CN-RCI")
        ax.set_title(f"{row.gene_name}\n({'masked' if row.gene_masked else 'visible'} at gene level)", fontsize=9)
    fig.suptitle("Robust compartment-switch exemplars (dots = all isoforms; large = tissue-dominant)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5.png")); plt.close(fig)


def fig6():
    # orthogonal validation: total RNA-seq vs polyA+ isoform CN-RCI (independent protocol)
    p = os.path.join(PROC, "orthogonal_total_rnaseq.tsv")
    if not os.path.exists(p):
        return
    d = pd.read_csv(p, sep="\t")
    fig, ax = plt.subplots(figsize=(4.2, 4.0))
    ax.scatter(d.rci_polyA, d.rci_total, s=10, alpha=0.5, color="#444")
    lim = [min(d.rci_polyA.min(), d.rci_total.min()) - 0.5, max(d.rci_polyA.max(), d.rci_total.max()) + 0.5]
    ax.plot(lim, lim, "--", color="#999", lw=1)
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("poly(A)+ isoform CN-RCI (15 cell lines)")
    ax.set_ylabel("total RNA-seq CN-RCI (K562+HepG2)")
    ax.set_title(f"Orthogonal protocol validation\n(n={KN['orthogonal_n_isoforms']}; Spearman ρ={KN['orthogonal_spearman_rho']}; "
                 f"sign agreement {KN['orthogonal_sign_agreement']})")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig6.png")); plt.close(fig)


if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    print("figures written:", sorted(os.listdir(FIG)))
