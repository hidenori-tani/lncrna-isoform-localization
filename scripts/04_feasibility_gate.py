"""Task 8: FEASIBILITY GATE.

Count switch events / genes where BOTH tissue-dominant isoforms have a localization
(non-NaN consensus CN-RCI). This is the GO/NO-GO number for the whole paper.

GO if genes_with_both_localized >= 50 (pre-registered in feasibility_report.md).
"""
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SE = os.path.join(ROOT, "data/raw/switch_exons.tsv")
LOC = os.path.join(ROOT, "data/processed/isoform_localization.tsv")
OUT = os.path.join(ROOT, "data/processed/switch_with_localization.tsv")


def strip_version(t):
    return str(t).split(".")[0]


def main():
    se = pd.read_csv(SE, sep="\t")
    loc = (
        pd.read_csv(LOC, sep="\t")
        .dropna(subset=["rci_consensus"])
        .set_index("transcript_id")
    )
    se["a"] = se.dom_isoform_A.map(strip_version)
    se["b"] = se.dom_isoform_B.map(strip_version)
    se["rci_a"] = se.a.map(loc.rci_consensus)
    se["rci_b"] = se.b.map(loc.rci_consensus)
    se["agree_a"] = se.a.map(loc.sign_agreement)
    se["agree_b"] = se.b.map(loc.sign_agreement)
    se["ncl_a"] = se.a.map(loc.n_cell_lines)
    se["ncl_b"] = se.b.map(loc.n_cell_lines)
    both = se.dropna(subset=["rci_a", "rci_b"]).copy()
    both.to_csv(OUT, sep="\t", index=False)

    print(f"switch_events_total={len(se)}")
    print(f"events_with_both_localized={len(both)}")
    print(f"genes_with_both_localized={both.gene_id.nunique()}")
    # 補助: dominant ペアが実際に異なる（=真のスイッチ）事象に限った数
    diff = both[both.a != both.b]
    print(f"events_with_both_localized_and_distinct_isoforms={len(diff)}")
    print(f"genes_with_distinct_localized_pairs={diff.gene_id.nunique()}")


if __name__ == "__main__":
    main()
