"""R3 validations (addresses 3-reviewer R2):
(a) INDEPENDENT QUANTIFIER: recompute isoform CN-RCI with kallisto (pseudo-alignment)
    on the SAME poly(A)+ ENCODE experiments and correlate with the RSEM (EM) consensus
    -> rules out shared-RSEM algorithmic bias.
(b) WITHIN-SAMPLE masking: define masking per cell line (a gene's isoforms occupy both
    compartments AND gene-level is neutral IN THE SAME sample), not via cross-cell-line
    median; report genes masked in >=1 cell line and the per-cell-line median.
(c) cross-cell-line CN-RCI variability (median per-isoform SD).
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd
from scipy import stats

from cswitch.localization import cn_rci, consensus_rci
from cswitch.keynumbers import set_number

QDIR = os.path.join(ROOT, "data/raw/encode_quant")
MAN = os.path.join(ROOT, "data/raw/encode_manifest.tsv")
SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")


def sv(x):
    return str(x).split(".")[0]


def fmt(f):
    with open(os.path.join(QDIR, f + ".tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
    if "transcript_id" in cols and "TPM" in cols:
        return "rsem"
    if "target_id" in cols and "tpm" in cols:
        return "kallisto"
    return "other"


def load_rsem(f, keep):
    df = pd.read_csv(os.path.join(QDIR, f + ".tsv"), sep="\t", usecols=["transcript_id", "gene_id", "TPM"])
    df = df[df.transcript_id.str.startswith("ENST")]
    df["tid"] = df.transcript_id.map(sv); df["gid"] = df.gene_id.map(sv)
    return df[df.gid.isin(keep)].groupby(["tid", "gid"], as_index=False)["TPM"].sum()


def load_kallisto(f, keep):
    df = pd.read_csv(os.path.join(QDIR, f + ".tsv"), sep="\t", usecols=["target_id", "tpm"])
    parts = df.target_id.str.split("|", expand=True)
    df["tid"] = parts[0].map(sv); df["gid"] = parts[1].map(sv)
    df = df[df.tid.str.startswith("ENST") & df.gid.isin(keep)]
    return df.rename(columns={"tpm": "TPM"}).groupby(["tid", "gid"], as_index=False)["TPM"].sum()


def per_cell_rci(man, keep, loader, fraction_col="fraction"):
    """returns long df transcript_id, gene_id, cell_line, rci (per cell line)."""
    rows = []
    for cell, g in man.groupby("cell_line"):
        nf = list(g[g[fraction_col] == "nucleus"].file)
        cf = list(g[g[fraction_col] == "cytoplasm"].file)
        if not nf or not cf:
            continue
        n = pd.concat([loader(f, keep) for f in nf]).groupby(["tid", "gid"])["TPM"].mean()
        c = pd.concat([loader(f, keep) for f in cf]).groupby(["tid", "gid"])["TPM"].mean()
        for key in n.index.intersection(c.index):
            tid, gid = key
            rows.append({"transcript_id": tid, "gene_id": gid, "cell_line": cell,
                         "rci": cn_rci(float(c[key]), float(n[key]), min_tpm=1.0)})
    return pd.DataFrame(rows)


def main():
    se = pd.read_csv(SE, sep="\t"); keep = set(se.gene_id.map(sv))
    man = pd.read_csv(MAN, sep="\t")
    man["f"] = man.file.map(fmt)
    rsem_man = man[man.f == "rsem"]; kal_man = man[man.f == "kallisto"]

    # ---- (a) kallisto vs RSEM (independent quantifier) ----
    rsem_long = per_cell_rci(rsem_man, keep, load_rsem)
    kal_long = per_cell_rci(kal_man, keep, load_kallisto)
    rsem_c = consensus_rci(rsem_long).rename(columns={"rci_consensus": "rci_rsem"})
    kal_c = consensus_rci(kal_long).rename(columns={"rci_consensus": "rci_kal"})
    m = rsem_c.dropna(subset=["rci_rsem"]).merge(
        kal_c.dropna(subset=["rci_kal"])[["transcript_id", "rci_kal"]], on="transcript_id")
    rho, _ = stats.spearmanr(m.rci_rsem, m.rci_kal)
    sign = float(np.mean(np.sign(m.rci_rsem) == np.sign(m.rci_kal)))
    set_number("quantifier_rsem_vs_kallisto_n", int(len(m)))
    set_number("quantifier_rsem_vs_kallisto_spearman", round(float(rho), 3))
    set_number("quantifier_rsem_vs_kallisto_sign_agreement", round(sign, 3))
    print(f"(a) RSEM vs kallisto (same polyA+ data, different algorithm): n={len(m)} "
          f"Spearman rho={rho:.3f} sign-agreement={sign:.3f}")

    # ---- (c) cross-cell-line variability (per-isoform SD across cell lines) ----
    sd = rsem_long.dropna(subset=["rci"]).groupby("transcript_id")["rci"].agg(["std", "count"])
    med_sd = float(sd[sd["count"] >= 2]["std"].median())
    set_number("cross_cellline_median_rci_sd", round(med_sd, 3))
    print(f"(c) cross-cell-line per-isoform CN-RCI median SD (>=2 lines): {med_sd:.3f}")

    # ---- (b) within-sample masking ----
    analysis_genes = set(pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"]).gene_id) & keep
    # per cell line: gene-level RCI = log2(sum cyto / sum nuc); need per-cell sums.
    # rebuild per-cell gene sums from RSEM
    gene_masked_cells = {}  # gid -> set of cells where masked
    gene_tested_cells = {}
    for cell, g in rsem_man.groupby("cell_line"):
        nf = list(g[g.fraction == "nucleus"].file); cf = list(g[g.fraction == "cytoplasm"].file)
        if not nf or not cf:
            continue
        nis = pd.concat([load_rsem(f, keep) for f in nf]).groupby(["tid", "gid"])["TPM"].mean()
        cis = pd.concat([load_rsem(f, keep) for f in cf]).groupby(["tid", "gid"])["TPM"].mean()
        # isoform RCI in this cell
        common = nis.index.intersection(cis.index)
        df = pd.DataFrame({"tid": [k[0] for k in common], "gid": [k[1] for k in common],
                           "rci": [cn_rci(float(cis[k]), float(nis[k]), min_tpm=1.0) for k in common]}).dropna()
        # gene-level RCI in this cell (sum over isoforms per gene per fraction)
        gn = nis.groupby(level="gid").sum(); gc = cis.groupby(level="gid").sum()
        for gid, grp in df.groupby("gid"):
            if gid not in analysis_genes or len(grp) < 2:
                continue
            gene_tested_cells.setdefault(gid, 0)
            gene_tested_cells[gid] += 1
            spans = (grp.rci.max() > 0) and (grp.rci.min() < 0)
            if gid in gn.index and gid in gc.index:
                grci = cn_rci(float(gc[gid]), float(gn[gid]), min_tpm=1.0)
                if spans and grci is not None and not np.isnan(grci) and abs(grci) < 1.0:
                    gene_masked_cells.setdefault(gid, set()).add(cell)
    genes_eval = [g for g in gene_tested_cells if gene_tested_cells[g] >= 1]
    masked_any = sum(1 for g in genes_eval if len(gene_masked_cells.get(g, set())) >= 1)
    # genes masked in >= half of the cell lines where they are testable
    masked_maj = sum(1 for g in genes_eval
                     if len(gene_masked_cells.get(g, set())) >= 0.5 * gene_tested_cells[g])
    set_number("within_sample_genes_evaluated", len(genes_eval))
    set_number("within_sample_masked_any_cellline", masked_any)
    set_number("within_sample_masked_any_frac", round(masked_any / len(genes_eval), 3))
    set_number("within_sample_masked_majority", masked_maj)
    set_number("within_sample_masked_majority_frac", round(masked_maj / len(genes_eval), 3))
    print(f"(b) within-sample masking: genes evaluated={len(genes_eval)} "
          f"masked in >=1 cell line={masked_any} ({masked_any/len(genes_eval)*100:.1f}%) "
          f"masked in majority={masked_maj} ({masked_maj/len(genes_eval)*100:.1f}%)")


if __name__ == "__main__":
    main()
