# scripts/01_setup_data.py
import os, shutil, sys, pandas as pd

SRC = "/Users/tanihidenori/claude-work/research/paper-lncrna-isoform-switching-dry"
DST = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FILES = {
    "data/processed/switch_exons.tsv":       "data/raw/switch_exons.tsv",
    "data/processed/visibility_longread.tsv":"data/raw/visibility_longread.tsv",
    "data/processed/lnc_isi_longread.tsv":   "data/raw/lnc_isi_longread.tsv",
    "data/raw/gencode.v26.lncRNA_transcripts.fa.gz": "data/raw/gencode.v26.lncRNA_transcripts.fa.gz",
    "data/raw/gencode.v26.long_noncoding_RNAs.gtf.gz": "data/raw/gencode.v26.long_noncoding_RNAs.gtf.gz",
}

def main():
    for src_rel, dst_rel in FILES.items():
        src = os.path.join(SRC, src_rel); dst = os.path.join(DST, dst_rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not os.path.exists(src):
            print(f"MISSING SOURCE: {src}"); sys.exit(1)
        shutil.copy2(src, dst)
        print(f"copied {dst}")
    se = pd.read_csv(os.path.join(DST, "data/raw/switch_exons.tsv"), sep="\t")
    need = {"gene_id","tissueA","tissueB","dom_isoform_A","dom_isoform_B"}
    assert need.issubset(se.columns), f"missing cols: {need - set(se.columns)}"
    print(f"switch_exons rows={len(se)} unique_genes={se.gene_id.nunique()} "
          f"unique_isoforms={pd.unique(se[['dom_isoform_A','dom_isoform_B']].values.ravel()).size}")

if __name__ == "__main__":
    main()
