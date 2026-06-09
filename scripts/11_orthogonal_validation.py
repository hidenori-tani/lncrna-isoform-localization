"""R2 orthogonal validation (addresses 3-reviewer R1: rule out polyA+ short-read artifact).

Independent library protocol: ENCODE TOTAL RNA-seq (rRNA-depleted; captures non-poly(A),
less 3' bias) subcellular fractionation for K562 and HepG2 — DISTINCT from the poly(A)+
RSEM data used to build the localization map. We recompute isoform CN-RCI from total
RNA-seq and correlate with the poly(A)+ consensus for switch-gene isoforms. High
concordance (and compartment-sign agreement) argues the localization calls are
protocol-robust, not a poly(A)+ quantification artifact.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd
import requests
from scipy import stats

from cswitch.localization import cn_rci, consensus_rci
from cswitch.keynumbers import set_number

BASE = "https://www.encodeproject.org"
HEAD = {"accept": "application/json"}
QDIR = os.path.join(ROOT, "data/raw/encode_total_quant")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
OUT = os.path.join(ROOT, "data/processed/orthogonal_total_rnaseq.tsv")
MIN_TPM = 1.0


def sv(x):
    return str(x).split(".")[0]


def search(fraction, cell):
    url = (f"{BASE}/search/?type=Experiment&status=released"
           f"&replicates.library.biosample.subcellular_fraction_term_name={fraction}"
           f"&assay_title=total+RNA-seq&biosample_ontology.term_name={cell}"
           f"&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
           f"&limit=all&format=json")
    r = requests.get(url, headers=HEAD, timeout=120); r.raise_for_status()
    return [e["accession"] for e in r.json().get("@graph", [])]


def quant_files(acc):
    r = requests.get(f"{BASE}/experiments/{acc}/?format=json", headers=HEAD, timeout=120)
    r.raise_for_status(); exp = r.json()
    out = []
    for f in exp.get("files", []):
        if (f.get("output_type") == "transcript quantifications"
                and f.get("file_format") == "tsv" and f.get("genome_annotation") == "V29"):
            out.append((f["accession"], BASE + f["href"]))
    return out


def download(acc, href):
    os.makedirs(QDIR, exist_ok=True)
    p = os.path.join(QDIR, acc + ".tsv")
    if not os.path.exists(p):
        r = requests.get(href, headers=HEAD, timeout=600); r.raise_for_status()
        open(p, "wb").write(r.content)
    return p


def is_rsem(p):
    with open(p) as fh:
        cols = fh.readline().split("\t")
    return "transcript_id" in cols and "TPM" in cols


def load(p):
    df = pd.read_csv(p, sep="\t", usecols=["transcript_id", "TPM"])
    df = df[df.transcript_id.str.startswith("ENST")]
    df["tid"] = df.transcript_id.map(sv)
    return df.groupby("tid", as_index=False)["TPM"].sum()


def main():
    rows = []
    for cell in ["K562", "HepG2"]:
        nuc, cyt = [], []
        for acc in search("nucleus", cell):
            for a, h in quant_files(acc):
                p = download(a, h)
                if is_rsem(p):
                    nuc.append(load(p))
        for acc in search("cytosol", cell):
            for a, h in quant_files(acc):
                p = download(a, h)
                if is_rsem(p):
                    cyt.append(load(p))
        if not nuc or not cyt:
            print(f"  {cell}: missing total RNA-seq fraction"); continue
        n = pd.concat(nuc).groupby("tid")["TPM"].mean()
        c = pd.concat(cyt).groupby("tid")["TPM"].mean()
        for tid in n.index.intersection(c.index):
            rows.append({"transcript_id": tid, "cell_line": cell,
                         "rci": cn_rci(float(c[tid]), float(n[tid]), min_tpm=MIN_TPM)})
    total = consensus_rci(pd.DataFrame(rows)).rename(columns={"rci_consensus": "rci_total"})

    poly = pd.read_csv(LOC, sep="\t").dropna(subset=["rci_consensus"])
    m = total.dropna(subset=["rci_total"]).merge(
        poly[["transcript_id", "rci_consensus"]], on="transcript_id").dropna()
    m = m.rename(columns={"rci_consensus": "rci_polyA"})
    m.to_csv(OUT, sep="\t", index=False)

    rho, prho = stats.spearmanr(m.rci_polyA, m.rci_total)
    r, pr = stats.pearsonr(m.rci_polyA, m.rci_total)
    sign_agree = float(np.mean(np.sign(m.rci_polyA) == np.sign(m.rci_total)))
    set_number("orthogonal_n_isoforms", int(len(m)))
    set_number("orthogonal_spearman_rho", round(float(rho), 3))
    set_number("orthogonal_pearson_r", round(float(r), 3))
    set_number("orthogonal_sign_agreement", round(sign_agree, 3))
    print(f"total-RNA-seq vs polyA+ CN-RCI: n={len(m)} isoforms (K562+HepG2)")
    print(f"  Spearman rho={rho:.3f} (p={prho:.2g}); Pearson r={r:.3f}; compartment-sign agreement={sign_agree:.3f}")


if __name__ == "__main__":
    main()
