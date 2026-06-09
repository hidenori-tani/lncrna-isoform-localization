"""Task 12: gene-level localization blind spot.

Central contribution after the honest pivot: gene-level localization (what LncATLAS-
style analysis reports = log2(sum cytoplasm TPM / sum nucleus TPM)) MASKS isoform-level
compartment heterogeneity. We recompute the true gene-level CN-RCI from ENCODE (summing
isoform TPMs per gene per fraction per cell line) and ask, for each switch gene:
does the gene-level value sit in the neutral zone while its isoforms span both compartments?
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd

from cswitch.localization import cn_rci, consensus_rci
from cswitch.blindspot import is_masked
from cswitch.keynumbers import set_number

QDIR = os.path.join(ROOT, "data/raw/encode_quant")
MAN = os.path.join(ROOT, "data/raw/encode_manifest.tsv")
SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
CALLS = os.path.join(ROOT, "data/processed/compartment_calls.tsv")
OUT = os.path.join(ROOT, "data/processed/blindspot.tsv")

MIN_TPM = 1.0
NEUTRAL = 1.0


def strip_version(x):
    return str(x).split(".")[0]


def is_rsem(file_acc):
    with open(os.path.join(QDIR, file_acc + ".tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
    return ("transcript_id" in cols) and ("TPM" in cols)


def load_rsem_gene(file_acc, keep_genes):
    """RSEM TSV -> 遺伝子(stripped)あたり TPM 合算（switch 遺伝子のみ）。"""
    df = pd.read_csv(os.path.join(QDIR, file_acc + ".tsv"),
                     sep="\t", usecols=["transcript_id", "gene_id", "TPM"])
    df = df[df["transcript_id"].str.startswith("ENST")]
    df["gid"] = df["gene_id"].map(strip_version)
    df = df[df["gid"].isin(keep_genes)]
    return df.groupby("gid", as_index=False)["TPM"].sum()


def gene_level_rci(man, keep_genes):
    rows = []
    for cell, g in man.groupby("cell_line"):
        nuc_files = list(g[g.fraction == "nucleus"].file)
        cyt_files = list(g[g.fraction == "cytoplasm"].file)
        if not nuc_files or not cyt_files:
            continue
        n = pd.concat([load_rsem_gene(f, keep_genes) for f in nuc_files]).groupby("gid")["TPM"].mean()
        c = pd.concat([load_rsem_gene(f, keep_genes) for f in cyt_files]).groupby("gid")["TPM"].mean()
        common = n.index.intersection(c.index)
        for gid in common:
            rows.append({"transcript_id": gid, "cell_line": cell,
                         "rci": cn_rci(float(c[gid]), float(n[gid]), min_tpm=MIN_TPM)})
    long = pd.DataFrame(rows)
    cons = consensus_rci(long).rename(columns={"transcript_id": "gene_id",
                                               "rci_consensus": "gene_rci"})
    return cons[["gene_id", "gene_rci"]]


def main():
    se = pd.read_csv(SE, sep="\t")
    keep_genes = set(se["gene_id"].map(strip_version))
    man = pd.read_csv(MAN, sep="\t")
    man = man[man["file"].map(is_rsem)].copy()

    gl = gene_level_rci(man, keep_genes).dropna(subset=["gene_rci"])
    iso = pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"])

    # 解析対象 = compartment_calls に出てくる遺伝子（両 dom 局在の163遺伝子）
    calls = pd.read_csv(CALLS, sep="\t")
    analysis_genes = set(calls.gene_id.map(strip_version))
    cs_genes = set(calls[calls.cs].gene_id.map(strip_version))

    recs = []
    for gid, grp in iso[iso.gene_id.isin(analysis_genes)].groupby("gene_id"):
        rcis = grp.rci_consensus.tolist()
        if len(rcis) < 2:
            continue
        grow = gl[gl.gene_id == gid]
        if grow.empty:
            continue
        g_rci = float(grow.gene_rci.iloc[0])
        spans_both = (max(rcis) > 0) and (min(rcis) < 0)
        recs.append({
            "gene_id": gid, "gene_rci": g_rci, "n_iso": len(rcis),
            "iso_rci_min": min(rcis), "iso_rci_max": max(rcis),
            "spans_both": spans_both,
            "masked": is_masked(g_rci, rcis, neutral=NEUTRAL),
            "is_compartment_switcher": gid in cs_genes,
        })
    bs = pd.DataFrame(recs)
    bs.to_csv(OUT, sep="\t", index=False)

    n_tested = len(bs)
    n_spans = int(bs.spans_both.sum())
    n_masked = int(bs.masked.sum())
    # compartment-switcher のうち gene-level が中立で隠れている割合
    csub = bs[bs.is_compartment_switcher]
    n_cs_masked = int(csub.masked.sum())
    set_number("blindspot_genes_tested", n_tested)
    set_number("blindspot_genes_span_both_compartments", n_spans)
    set_number("blindspot_genes_masked", n_masked)
    set_number("blindspot_masked_fraction", float(n_masked / n_tested) if n_tested else 0.0)
    set_number("blindspot_cs_genes_masked", n_cs_masked)
    set_number("blindspot_cs_genes_total", int(len(csub)))

    print(f"genes tested (>=2 localized isoforms + gene-level): {n_tested}")
    print(f"  isoforms span both compartments: {n_spans} ({n_spans/n_tested*100:.1f}%)")
    print(f"  MASKED (gene-level |RCI|<{NEUTRAL} but isoforms span both): {n_masked} ({n_masked/n_tested*100:.1f}%)")
    print(f"  among compartment-switchers ({len(csub)}): masked={n_cs_masked} ({n_cs_masked/max(len(csub),1)*100:.1f}%)")


if __name__ == "__main__":
    main()
