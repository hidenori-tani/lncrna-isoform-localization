"""Task 14: Alu/SINE mechanism.

Even though tissue switching is compartment-agnostic, the localization DIFFERENCE
between two isoforms should be sequence-encoded if it is real (not noise). Test:
do isoform pairs that occupy DIFFERENT compartments (compartment-switch events) differ
more in Alu content than concordant pairs, controlling for transcript-length difference?
A positive result argues the compartment differences are biologically real (Alu drives
nuclear retention), strengthening the blind-spot thesis against a 'pure noise' critique.
"""
import gzip
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import numpy as np
import pandas as pd

from cswitch.mechanism import alu_content_diff, length_matched_compare
from cswitch.keynumbers import set_number

GTF = os.path.join(ROOT, "data/raw/gencode.v26.long_noncoding_RNAs.gtf.gz")
RMSK = os.path.join(ROOT, "data/raw/rmsk.txt.gz")
CALLS = os.path.join(ROOT, "data/processed/compartment_calls.tsv")
OUT_ISO = os.path.join(ROOT, "data/processed/isoform_alu.tsv")
OUT_PAIR = os.path.join(ROOT, "data/processed/mechanism_pairs.tsv")


def strip_version(x):
    return str(x).split(".")[0]


def parse_exons(needed):
    """needed: set of stripped transcript ids -> dict tid -> list[(chr, s0, e0)] (0-based half-open)."""
    exons = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            attrs = f[8]
            i = attrs.find('transcript_id "')
            if i < 0:
                continue
            tid = strip_version(attrs[i + 15:attrs.find('"', i + 15)])
            if tid not in needed:
                continue
            exons.setdefault(tid, []).append((f[0], int(f[3]) - 1, int(f[4])))
    return exons


def load_alu():
    """rmsk -> dict chr -> (starts[sorted], ends[reordered]) for repFamily=='Alu'."""
    df = pd.read_csv(RMSK, sep="\t", header=None,
                     usecols=[5, 6, 7, 12], names=["chr", "start", "end", "fam"],
                     compression="gzip")
    df = df[df.fam == "Alu"]
    out = {}
    for chrom, g in df.groupby("chr"):
        g = g.sort_values("start")
        out[chrom] = (g.start.to_numpy(), g.end.to_numpy())
    return out


def alu_bp_for_exon(alu, chrom, s, e):
    if chrom not in alu:
        return 0
    starts, ends = alu[chrom]
    # candidate Alu intervals with start < e (binary search upper bound)
    hi = np.searchsorted(starts, e, side="left")
    if hi == 0:
        return 0
    cs = starts[:hi]
    ce = ends[:hi]
    # overlap bp with [s, e): sum max(0, min(e, ce) - max(s, cs)), only where ce > s
    ov = np.minimum(e, ce) - np.maximum(s, cs)
    return int(ov[ov > 0].sum())


def main():
    calls = pd.read_csv(CALLS, sep="\t")
    needed = set(calls.a) | set(calls.b)
    print(f"needed transcripts: {len(needed)}")

    exons = parse_exons(needed)
    print(f"transcripts with exon annotation: {len(exons)}")
    alu = load_alu()

    rows = []
    for tid, exl in exons.items():
        length = sum(e - s for _, s, e in exl)
        abp = sum(alu_bp_for_exon(alu, c, s, e) for c, s, e in exl)
        rows.append({"transcript_id": tid, "length": length,
                     "alu_bp": abp, "alu_frac": (abp / length if length else np.nan)})
    iso = pd.DataFrame(rows)
    iso.to_csv(OUT_ISO, sep="\t", index=False)
    af = iso.set_index("transcript_id")["alu_frac"]
    ln = iso.set_index("transcript_id")["length"]

    calls = calls.copy()
    calls["alu_a"] = calls.a.map(af); calls["alu_b"] = calls.b.map(af)
    calls["len_a"] = calls.a.map(ln); calls["len_b"] = calls.b.map(ln)
    pair = calls.dropna(subset=["alu_a", "alu_b", "len_a", "len_b"]).copy()
    pair["alu_diff"] = [alu_content_diff(a, b) for a, b in zip(pair.alu_a, pair.alu_b)]
    pair["len_diff"] = (pair.len_a - pair.len_b).abs()
    pair["group"] = np.where(pair.cs, "cs", "concordant")
    pair.to_csv(OUT_PAIR, sep="\t", index=False)

    res = length_matched_compare(pair, value="alu_diff", group="group", covariate="len_diff")
    set_number("mechanism_alu_diff_p", res["p_value"])
    set_number("mechanism_alu_diff_median_cs", res["median_cs"])
    set_number("mechanism_alu_diff_median_concordant", res["median_concordant"])
    set_number("mechanism_n_pairs", int(len(pair)))
    set_number("mechanism_n_cs_pairs", int((pair.group == "cs").sum()))

    print(f"pairs with Alu+length: {len(pair)} (cs={int((pair.group=='cs').sum())})")
    print(f"Alu diff (cs vs concordant): median_cs={res['median_cs']:.4f} "
          f"median_concordant={res['median_concordant']:.4f} p(one-sided)={res['p_value']:.4g}")


if __name__ == "__main__":
    main()
