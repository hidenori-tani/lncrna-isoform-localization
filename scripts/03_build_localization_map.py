"""Task 7: build per-isoform subcellular localization map (CN-RCI) from ENCODE
fractionation RSEM transcript quantifications.

Design decisions (data-driven, 2026-06-09):
- Downloaded "transcript quantifications / V29" files come in TWO formats: RSEM
  (transcript_id, gene_id, TPM) and kallisto (target_id, tpm). RSEM exists for ALL
  15 cell lines in BOTH fractions, so we standardize on RSEM only (avoid tool-mixing
  bias in the cross-cell-line consensus).
- RSEM ENST rows carry ENST.version (transcript) and ENSG.version (gene). We strip
  versions to bridge ENCODE V29 -> reused GENCODE v26 switch calls (stable ID match).
- We restrict to transcripts of the 268 switch genes (from switch_exons.tsv): this is
  exactly the universe needed downstream (feasibility gate, compartment calls,
  within-gene permutation null, gene-level blind spot) and keeps it fast.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)  # allow `import cswitch` when run as `python scripts/...`

import pandas as pd

from cswitch.localization import cn_rci, consensus_rci

QDIR = os.path.join(ROOT, "data/raw/encode_quant")
MAN = os.path.join(ROOT, "data/raw/encode_manifest.tsv")
SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
OUT = os.path.join(ROOT, "data/processed/isoform_localization.tsv")

MIN_TPM = 1.0  # localization は MIN_TPM 以上の検出を要求（低発現ノイズ抑制・cn_rci 既定0.1より厳格）


def strip_version(x):
    return str(x).split(".")[0]


def is_rsem(file_acc):
    with open(os.path.join(QDIR, file_acc + ".tsv")) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
    return ("transcript_id" in cols) and ("TPM" in cols)


def load_rsem(file_acc, keep_genes):
    """RSEM TSV を読み、switch 遺伝子に属する ENST 行のみ返す。
    返り値: DataFrame[tid, gid, TPM]（tid/gid は version 落とし）。"""
    df = pd.read_csv(
        os.path.join(QDIR, file_acc + ".tsv"),
        sep="\t",
        usecols=["transcript_id", "gene_id", "TPM"],
    )
    df = df[df["transcript_id"].str.startswith("ENST")]
    df["tid"] = df["transcript_id"].map(strip_version)
    df["gid"] = df["gene_id"].map(strip_version)
    df = df[df["gid"].isin(keep_genes)]
    # 同一 stable id が複数 version で現れたら合算
    return df.groupby(["tid", "gid"], as_index=False)["TPM"].sum()


def fraction_tpm(files, keep_genes):
    """複数 RSEM ファイル（同一細胞株・同一分画の複製）を平均した tid->TPM。"""
    parts = [load_rsem(f, keep_genes) for f in files]
    cat = pd.concat(parts, ignore_index=True)
    # 複製間は平均（gid も保持）
    g = cat.groupby(["tid", "gid"], as_index=False)["TPM"].mean()
    return g


def main():
    se = pd.read_csv(SE, sep="\t")
    keep_genes = set(se["gene_id"].map(strip_version))
    print(f"switch genes (stripped): {len(keep_genes)}")

    man = pd.read_csv(MAN, sep="\t")
    man = man[man["file"].map(is_rsem)].copy()
    print(f"RSEM files: {len(man)}")

    rows = []
    tid2gid = {}
    for cell, g in man.groupby("cell_line"):
        nuc_files = list(g[g.fraction == "nucleus"].file)
        cyt_files = list(g[g.fraction == "cytoplasm"].file)
        if not nuc_files or not cyt_files:
            continue
        n = fraction_tpm(nuc_files, keep_genes).set_index("tid")
        c = fraction_tpm(cyt_files, keep_genes).set_index("tid")
        common = n.index.intersection(c.index)
        for tid in common:
            r = cn_rci(float(c.loc[tid, "TPM"]), float(n.loc[tid, "TPM"]), min_tpm=MIN_TPM)
            rows.append({"transcript_id": tid, "cell_line": cell, "rci": r})
            tid2gid[tid] = n.loc[tid, "gid"]
    long = pd.DataFrame(rows)
    print(f"per-(transcript,cell_line) RCI rows: {len(long)}")

    cons = consensus_rci(long)
    cons["gene_id"] = cons["transcript_id"].map(tid2gid)
    cons = cons[["transcript_id", "gene_id", "rci_consensus", "n_cell_lines", "sign_agreement"]]
    cons.to_csv(OUT, sep="\t", index=False)

    nloc = int(cons.rci_consensus.notna().sum())
    n2 = int((cons.n_cell_lines >= 2).sum())
    print(f"isoforms with localization: {nloc}  (>=2 cell lines: {n2})")
    print(f"unique switch genes covered: {cons.loc[cons.rci_consensus.notna(),'gene_id'].nunique()}")


if __name__ == "__main__":
    main()
