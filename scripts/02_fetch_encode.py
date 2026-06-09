# scripts/02_fetch_encode.py
import requests, json, os, sys, csv

BASE = "https://www.encodeproject.org"
HEAD = {"accept": "application/json"}
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data/raw")

def search_experiments(fraction):
    url = (f"{BASE}/search/?type=Experiment&status=released"
           f"&replicates.library.biosample.subcellular_fraction_term_name={fraction}"
           f"&assay_title=polyA+plus+RNA-seq&assay_title=total+RNA-seq"
           f"&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens"
           f"&limit=all&format=json")
    r = requests.get(url, headers=HEAD, timeout=120); r.raise_for_status()
    return r.json().get("@graph", [])

def transcript_quant_files(acc):
    r = requests.get(f"{BASE}/experiments/{acc}/?format=json", headers=HEAD, timeout=120)
    r.raise_for_status(); exp = r.json()
    rows = []
    for f in exp.get("files", []):
        # 単一 annotation (V29) に固定。固定しないと V19/V24/V29 が混ざり同一 isoform が重複する。
        if (f.get("output_type") == "transcript quantifications"
                and f.get("file_format") == "tsv"
                and f.get("genome_annotation") == "V29"):
            rows.append({
                "experiment": acc,
                "biosample": exp.get("biosample_summary",""),
                "cell_line": exp.get("biosample_ontology",{}).get("term_name",""),
                "file": f["accession"],
                "href": BASE + f["href"],
                "assembly": f.get("assembly",""),
                "genome_annotation": f.get("genome_annotation",""),
            })
    return rows

def main():
    manifest = []
    # ENCODE の細胞質画分ラベルは "cytosol"。cytosol を内部ラベル "cytoplasm" に正規化。
    for frac in ["nucleus", "cytosol"]:
        for exp in search_experiments(frac):
            acc = exp["accession"]
            for row in transcript_quant_files(acc):
                row["fraction"] = "nucleus" if frac == "nucleus" else "cytoplasm"
                manifest.append(row)
    if not manifest:
        print("NO FRACTION EXPERIMENTS FOUND — schema を docs/encode_manifest.md に手動確認"); sys.exit(2)
    path = os.path.join(OUT, "encode_manifest.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(manifest)
    cells = sorted({m["cell_line"] for m in manifest})
    print(f"manifest rows={len(manifest)} cell_lines={len(cells)}: {cells}")

def download(manifest_path, dst_dir, limit=None):
    import pandas as pd
    os.makedirs(dst_dir, exist_ok=True)
    man = pd.read_csv(manifest_path, sep="\t")
    if limit is not None:
        man = man.head(limit)
    for _, row in man.iterrows():
        out = os.path.join(dst_dir, f"{row['file']}.tsv")
        if os.path.exists(out):
            continue
        r = requests.get(row["href"], headers=HEAD, timeout=600); r.raise_for_status()
        with open(out, "wb") as fh:
            fh.write(r.content)
        print(f"downloaded {out} ({len(r.content)} bytes)")

if __name__ == "__main__":
    if "--download" in sys.argv:
        ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        lim = None
        if "--limit" in sys.argv:
            lim = int(sys.argv[sys.argv.index("--limit") + 1])
        download(os.path.join(ROOT, "data/raw/encode_manifest.tsv"),
                 os.path.join(ROOT, "data/raw/encode_quant"), limit=lim)
    else:
        main()
