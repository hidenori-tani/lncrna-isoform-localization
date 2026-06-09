"""Task 15: case studies — robust exemplars of isoform-resolved compartment divergence
that gene-level analysis masks. Selection is honest: require BOTH isoforms robustly
localized (sign_agreement>=0.8, n_cell_lines>=3), large |dRCI|, and gene-level masked.
Gene symbols from GENCODE GTF. Deep literature/PMID verification is deferred to QA.
"""
import gzip
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd

GTF = os.path.join(ROOT, "data/raw/gencode.v26.long_noncoding_RNAs.gtf.gz")
CALLS = os.path.join(ROOT, "data/processed/compartment_calls.tsv")
BS = os.path.join(ROOT, "data/processed/blindspot.tsv")
OUT = os.path.join(ROOT, "data/processed/case_studies.tsv")


def strip_version(x):
    return str(x).split(".")[0]


def gene_symbols():
    sym = {}
    with gzip.open(GTF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            a = f[8]
            gi = a.find('gene_id "'); gn = a.find('gene_name "')
            if gi < 0 or gn < 0:
                continue
            gid = strip_version(a[gi + 9:a.find('"', gi + 9)])
            name = a[gn + 11:a.find('"', gn + 11)]
            sym[gid] = name
    return sym


def main():
    calls = pd.read_csv(CALLS, sep="\t")
    bs = pd.read_csv(BS, sep="\t")
    masked_genes = set(bs[bs.masked].gene_id)
    sym = gene_symbols()

    c = calls[calls.cs].copy()
    c["gid"] = c.gene_id.map(strip_version)
    c = c[(c.agree_a >= 0.8) & (c.agree_b >= 0.8) & (c.ncl_a >= 3) & (c.ncl_b >= 3)]
    c["abs_drci"] = c.delta_rci.abs()
    c = c[c.abs_drci >= 2.0]
    c["gene_masked"] = c.gid.isin(masked_genes)
    c["gene_name"] = c.gid.map(sym)
    c = c.sort_values("abs_drci", ascending=False)
    top = c.drop_duplicates("gid").head(9)[
        ["gid", "gene_name", "tissueA", "tissueB", "dom_isoform_A", "dom_isoform_B",
         "rci_a", "rci_b", "delta_rci", "agree_a", "agree_b", "ncl_a", "ncl_b", "gene_masked"]
    ]
    top.to_csv(OUT, sep="\t", index=False)
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print(f"robust compartment-switch genes (|dRCI|>=2, both agree>=0.8, ncl>=3): {c.gid.nunique()}")
    print("\nTop candidates (gene-level masked flagged):")
    print(top[["gid", "gene_name", "tissueA", "tissueB", "rci_a", "rci_b", "delta_rci", "gene_masked"]].to_string(index=False))


if __name__ == "__main__":
    main()
