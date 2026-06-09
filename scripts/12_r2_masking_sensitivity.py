"""R2: (a) abundance-aware masking (dominant isoforms only) and
       (b) min_TPM sensitivity of the 'span-both-compartments' fraction.
Addresses reviewer asks: masking should be abundance-aware; CN-RCI sensitivity.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd

from cswitch.localization import cn_rci, consensus_rci
from cswitch.keynumbers import set_number

QDIR = os.path.join(ROOT, "data/raw/encode_quant")
MAN = os.path.join(ROOT, "data/raw/encode_manifest.tsv")
SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
BS = os.path.join(ROOT, "data/processed/blindspot.tsv")


def sv(x):
    return str(x).split(".")[0]


# ---------- (a) abundance-aware masking ----------
def abundance_aware_masking():
    se = pd.read_csv(SE, sep="\t"); se["gid"] = se.gene_id.map(sv)
    loc = pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"]).set_index("transcript_id")["rci_consensus"]
    bs = pd.read_csv(BS, sep="\t").set_index("gene_id")  # has gene_rci per gene (stripped)
    n_tested = n_span_dom = n_masked_dom = 0
    for gid, g in se.groupby("gid"):
        if gid not in bs.index:
            continue
        dom = set(g.dom_isoform_A.map(sv)) | set(g.dom_isoform_B.map(sv))
        rcis = [float(loc[t]) for t in dom if t in loc.index]
        if len(rcis) < 2:
            continue
        n_tested += 1
        spans = (max(rcis) > 0) and (min(rcis) < 0)
        if spans:
            n_span_dom += 1
            if abs(float(bs.loc[gid, "gene_rci"])) < 1.0:
                n_masked_dom += 1
    set_number("blindspot_dominant_genes_tested", n_tested)
    set_number("blindspot_dominant_span_both", n_span_dom)
    set_number("blindspot_dominant_span_both_frac", round(n_span_dom / n_tested, 3))
    set_number("blindspot_dominant_masked", n_masked_dom)
    set_number("blindspot_dominant_masked_frac", round(n_masked_dom / n_tested, 3))
    print(f"[abundance-aware] genes with >=2 localized DOMINANT isoforms: {n_tested}")
    print(f"  dominant isoforms span both compartments: {n_span_dom} ({n_span_dom/n_tested*100:.1f}%)")
    print(f"  masked (gene-level neutral + dominant span both): {n_masked_dom} ({n_masked_dom/n_tested*100:.1f}%)")


# ---------- (b) min_TPM sensitivity ----------
def is_rsem(f):
    with open(os.path.join(QDIR, f + ".tsv")) as fh:
        cols = fh.readline().split("\t")
    return "transcript_id" in cols and "TPM" in cols


def load(f, keep):
    df = pd.read_csv(os.path.join(QDIR, f + ".tsv"), sep="\t", usecols=["transcript_id", "gene_id", "TPM"])
    df = df[df.transcript_id.str.startswith("ENST")]
    df["tid"] = df.transcript_id.map(sv); df["gid"] = df.gene_id.map(sv)
    df = df[df.gid.isin(keep)]
    return df.groupby(["tid", "gid"], as_index=False)["TPM"].sum()


def min_tpm_sensitivity():
    se = pd.read_csv(SE, sep="\t"); keep = set(se.gene_id.map(sv))
    man = pd.read_csv(MAN, sep="\t"); man = man[man.file.map(is_rsem)]
    # read raw nuc/cyt TPM per cell line ONCE (store), then apply floors
    per_cell = {}
    for cell, g in man.groupby("cell_line"):
        nf = list(g[g.fraction == "nucleus"].file); cf = list(g[g.fraction == "cytoplasm"].file)
        if not nf or not cf:
            continue
        n = pd.concat([load(f, keep) for f in nf]).groupby(["tid", "gid"])["TPM"].mean()
        c = pd.concat([load(f, keep) for f in cf]).groupby(["tid", "gid"])["TPM"].mean()
        per_cell[cell] = (n, c)
    out = {}
    for floor in [0.5, 1.0, 2.0]:
        rows = []
        for cell, (n, c) in per_cell.items():
            common = n.index.intersection(c.index)
            for (tid, gid) in common:
                rows.append({"transcript_id": tid, "gid": gid, "cell_line": cell,
                             "rci": cn_rci(float(c[(tid, gid)]), float(n[(tid, gid)]), min_tpm=floor)})
        long = pd.DataFrame(rows)
        cons = consensus_rci(long.rename(columns={"transcript_id": "transcript_id"})).dropna(subset=["rci_consensus"])
        # attach gid
        tid2gid = long.drop_duplicates("transcript_id").set_index("transcript_id")["gid"]
        cons["gid"] = cons.transcript_id.map(tid2gid)
        n_t = n_s = 0
        for gid, gg in cons.groupby("gid"):
            r = gg.rci_consensus.tolist()
            if len(r) < 2:
                continue
            n_t += 1
            if max(r) > 0 and min(r) < 0:
                n_s += 1
        out[str(floor)] = round(n_s / n_t, 3) if n_t else None
        print(f"  min_TPM={floor}: genes={n_t} span_both={n_s} ({(n_s/n_t*100 if n_t else 0):.1f}%)")
    set_number("blindspot_span_both_by_min_tpm", out)


if __name__ == "__main__":
    print("=== (a) abundance-aware masking ===")
    abundance_aware_masking()
    print("=== (b) min_TPM sensitivity (span-both fraction) ===")
    min_tpm_sensitivity()
