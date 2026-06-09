"""R2 revision analyses (addressing 3-reviewer R1 major points):
- mapping fidelity: how many v26 dominant isoforms / genes have ENCODE V29 localization
- expression/dominance-matched null: null pool = DOMINANT isoforms only (well-expressed,
  not minor pre-mRNA noise) instead of all localized isoforms
- representativeness: retained 163 vs input 268 genes (n dominant isoforms, lnc-ISI)
All numbers -> key_numbers.json.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd

from cswitch.compartment_switching import permutation_null_gene_fraction
from cswitch.keynumbers import set_number

SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
SWL = os.path.join(ROOT, "data/processed/switch_with_localization.tsv")
CALLS = os.path.join(ROOT, "data/processed/compartment_calls.tsv")
LNCISI = os.path.join(ROOT, "data/raw/lnc_isi_longread.tsv")


def sv(x):
    return str(x).split(".")[0]


def main():
    se = pd.read_csv(SE, sep="\t")
    loc = pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"])
    localized = set(loc.transcript_id)

    # ---- mapping fidelity ----
    dom_iso = set(se.dom_isoform_A.map(sv)) | set(se.dom_isoform_B.map(sv))
    dom_iso_localized = dom_iso & localized
    map_rate = len(dom_iso_localized) / len(dom_iso)
    set_number("mapping_dominant_isoforms_total", len(dom_iso))
    set_number("mapping_dominant_isoforms_localized", len(dom_iso_localized))
    set_number("mapping_dominant_isoform_rate", round(map_rate, 4))
    print(f"dominant isoforms: {len(dom_iso)}  localized(V29-bridged): {len(dom_iso_localized)} "
          f"({map_rate*100:.1f}%)")

    # ---- dominance-matched null (pool = DOMINANT isoforms only) ----
    calls = pd.read_csv(CALLS, sep="\t")
    d = calls.copy()
    d["gid"] = d.gene_id.map(sv)
    gene_events = d.groupby("gid").size().to_dict()
    frac = float(d.groupby("gid")["cs"].any().mean())  # observed fraction (== compartment_switcher_gene_fraction)
    # pool per gene = that gene's DOMINANT isoforms that are localized
    se2 = se.copy(); se2["gid"] = se2.gene_id.map(sv)
    loc_rci = loc.set_index("transcript_id")["rci_consensus"]
    gene_pool = {}
    for gid, g in se2.groupby("gid"):
        iso = set(g.dom_isoform_A.map(sv)) | set(g.dom_isoform_B.map(sv))
        rcis = [float(loc_rci[t]) for t in iso if t in loc_rci.index]
        if len(rcis) >= 2:
            gene_pool[gid] = np.array(rcis)
    _, p_dom, null_dom = permutation_null_gene_fraction(
        frac, gene_events, gene_pool, threshold=1.0, B=10000, seed=42)
    set_number("compartment_null_dominant_pool_mean", float(null_dom.mean()))
    set_number("compartment_null_dominant_pool_p", p_dom)
    set_number("compartment_null_dominant_pool_genes", len(gene_pool))
    print(f"dominance-matched null: observed={frac:.3f} null_mean={null_dom.mean():.3f} "
          f"p={p_dom:.4g} (pool genes={len(gene_pool)})")

    # ---- representativeness: retained 163 vs input 268 ----
    swl = pd.read_csv(SWL, sep="\t")
    retained = set(swl.gene_id.map(sv))
    allg = set(se.gene_id.map(sv))
    # n dominant isoforms per gene
    nd = se2.groupby("gid").apply(
        lambda g: len(set(g.dom_isoform_A) | set(g.dom_isoform_B)), include_groups=False)
    lncisi = pd.read_csv(LNCISI, sep="\t"); lncisi["gid"] = lncisi.gene_id.map(sv)
    isi = lncisi.set_index("gid")["lnc_isi"]
    def summ(genes, name):
        nds = [nd[g] for g in genes if g in nd.index]
        isis = [isi[g] for g in genes if g in isi.index]
        print(f"  {name} (n={len(genes)}): median dom-isoforms={np.median(nds):.1f}, "
              f"median lnc-ISI={np.median(isis):.3f}")
        return float(np.median(nds)), float(np.median(isis))
    rn, ri = summ(retained, "retained-163")
    an, ai = summ(allg, "input-268")
    set_number("representativeness_retained_median_dom_isoforms", rn)
    set_number("representativeness_input_median_dom_isoforms", an)
    set_number("representativeness_retained_median_lnc_isi", ri)
    set_number("representativeness_input_median_lnc_isi", ai)


if __name__ == "__main__":
    main()
