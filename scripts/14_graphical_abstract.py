"""Graphical abstract (CBC requires one at submission).

Elsevier spec: image 531 x 1328 px (h x w) or proportionally more, readable at
5 x 13 cm. We render at 2x (1062 x 2656 px) as PNG + PDF.

All numbers derive from key_numbers.json (no hardcoded results).
Concept: the SAME switching lncRNA gene looks neutral at gene level (masked)
but its isoforms occupy opposite compartments at isoform resolution.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

FIG = os.path.join(ROOT, "figures")
KN = json.load(open(os.path.join(ROOT, "key_numbers.json")))

span_both_pct = round(KN["blindspot_genes_span_both_compartments"] / KN["blindspot_genes_tested"] * 100)
masked_pct = round(KN["blindspot_masked_fraction"] * 100)
rho = KN["orthogonal_spearman_rho"]

NUC = "#2c6fbb"   # nuclear
CYT = "#d1495b"   # cytoplasmic
GREY = "#8a8a8a"
INK = "#222222"

# 1328 x 531 (w x h) aspect = 2.5:1; render at 2x via dpi
fig = plt.figure(figsize=(13.28, 5.31), dpi=200)
fig.patch.set_facecolor("white")
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis("off")

# ---- headline ----
ax.text(50, 37.5, "lncRNA subcellular localization is an isoform-resolved property",
        ha="center", va="center", fontsize=18, fontweight="bold", color=INK)
ax.text(50, 33.8, "masked by gene-level analysis", ha="center", va="center",
        fontsize=18, fontweight="bold", color=INK)

# ===== LEFT: gene-level view (masked) =====
lx = 19
ax.text(lx, 29, "Gene-level view", ha="center", fontsize=14, fontweight="bold", color=INK)
# axis
ax.annotate("", xy=(lx + 11, 19), xytext=(lx - 11, 19),
            arrowprops=dict(arrowstyle="-", color=INK, lw=1.5))
ax.text(lx - 11, 16.7, "nuclear", ha="center", fontsize=10, color=NUC)
ax.text(lx + 11, 16.7, "cytoplasmic", ha="center", fontsize=10, color=CYT)
ax.axvline  # noop
ax.plot([lx, lx], [18.3, 19.7], color=INK, lw=1, ls="--")  # zero tick
# single neutral dot = gene aggregate
ax.scatter([lx], [19], s=900, color=GREY, edgecolor=INK, zorder=5)
ax.text(lx, 19, "gene", ha="center", va="center", fontsize=9, color="white", fontweight="bold")
ax.text(lx, 13.5, "one neutral value\n→ looks unlocalized", ha="center", fontsize=11, color=INK)
ax.text(lx, 7.5, f"masked in {masked_pct}% of genes", ha="center", fontsize=13,
        fontweight="bold", color="#b03050")

# ===== arrow =====
arr = FancyArrowPatch((36, 19), (47, 19), arrowstyle="-|>", mutation_scale=28,
                      color=INK, lw=2.2)
ax.add_patch(arr)
ax.text(41.5, 22, "resolve\nisoforms", ha="center", fontsize=10, color=INK)

# ===== RIGHT: isoform-resolved view (split) =====
rx = 73
ax.text(rx, 29, "Isoform-resolved view", ha="center", fontsize=14, fontweight="bold", color=INK)
ax.annotate("", xy=(rx + 18, 19), xytext=(rx - 18, 19),
            arrowprops=dict(arrowstyle="-", color=INK, lw=1.5))
ax.text(rx - 18, 16.7, "nuclear", ha="center", fontsize=10, color=NUC)
ax.text(rx + 18, 16.7, "cytoplasmic", ha="center", fontsize=10, color=CYT)
ax.plot([rx, rx], [18.3, 19.7], color=INK, lw=1, ls="--")
# isoforms split to both compartments
rng = np.random.default_rng(0)
nuc_x = rx - (6 + rng.random(4) * 9)
cyt_x = rx + (6 + rng.random(4) * 9)
ax.scatter(nuc_x, 19 + rng.normal(0, 0.0, 4), s=320, color=NUC, edgecolor=INK, zorder=5)
ax.scatter(cyt_x, 19 + rng.normal(0, 0.0, 4), s=320, color=CYT, edgecolor=INK, zorder=5)
ax.text(rx, 13.5, "isoforms occupy\nopposite compartments", ha="center", fontsize=11, color=INK)
ax.text(rx, 7.5, f"{span_both_pct}% of switching lncRNA genes", ha="center", fontsize=13,
        fontweight="bold", color="#1f5fae")

# ---- footer: validation ----
ax.text(50, 2.4, f"Measured from ENCODE subcellular fractionation (15 cell lines); "
                 f"reproduced by an independent total RNA-seq protocol (Spearman ρ = {rho})",
        ha="center", fontsize=9.5, color=GREY)

png = os.path.join(FIG, "graphical_abstract.png")
pdf = os.path.join(FIG, "graphical_abstract.pdf")
fig.savefig(png, dpi=200, facecolor="white")
fig.savefig(pdf, facecolor="white")
plt.close(fig)

# verify pixel size
try:
    from PIL import Image
    w, h = Image.open(png).size
    print(f"graphical abstract: {png}  {w}x{h}px (need >= 1328x531 w x h)")
    assert w >= 1328 and h >= 531, "graphical abstract too small"
except ImportError:
    print(f"graphical abstract written: {png} (install Pillow to verify size)")
