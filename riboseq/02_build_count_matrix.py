#!/usr/bin/env python3
"""Combine per-sample featureCounts CDS tables of one batch into a gene x sample matrix.

Usage: riboseq/02_build_count_matrix.py --batch <name> --indir <processed_dir> --outdir <dir>
Input:  <processed_dir>/featureCounts_cds_<sample>.dedup.tsv
Output: <dir>/counts_<batch>.tsv, <dir>/coldata_<batch>.tsv, <dir>/gene_info.tsv
"""
import argparse
import re
from pathlib import Path

SAMPLE_RE = re.compile(r"(?:(?P<prefix>[AB]\d+)_)?(?P<fraction>80S|disome)_(?P<condition>[A-Za-z0-9]+?)R(?P<replicate>\d)_(?P<barcode>[ACGT]+)$")


def load_counts(tsv):
    info, counts = {}, {}
    with open(tsv) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            if header is None:
                header = line.rstrip("\n").split("\t")
                continue
            p = line.rstrip("\n").split("\t")
            info[p[0]] = p[1:6]
            counts[p[0]] = int(p[6])
    return info, counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", required=True)
    ap.add_argument("--indir", required=True, type=Path)
    ap.add_argument("--outdir", required=True, type=Path)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    samples, gene_info, counts_by_sample = [], None, {}
    for tsv in sorted(a.indir.glob("featureCounts_cds_*.dedup.tsv")):
        sample = tsv.name[len("featureCounts_cds_"):-len(".dedup.tsv")]
        info, counts = load_counts(tsv)
        gene_info = gene_info or info
        samples.append(sample)
        counts_by_sample[sample] = counts
    genes = list(gene_info)

    with open(a.outdir / f"counts_{a.batch}.tsv", "w") as out:
        out.write("gene_name\t" + "\t".join(samples) + "\n")
        for g in genes:
            out.write(g + "\t" + "\t".join(str(counts_by_sample[s].get(g, 0)) for s in samples) + "\n")

    with open(a.outdir / f"coldata_{a.batch}.tsv", "w") as out:
        out.write("sample\tprefix\tfraction\tcondition\treplicate\tbarcode\n")
        for s in samples:
            m = SAMPLE_RE.match(s)
            d = m.groupdict() if m else {}
            out.write("\t".join([s] + [d.get(k) or "NA" for k in ("prefix", "fraction", "condition", "replicate", "barcode")]) + "\n")

    with open(a.outdir / "gene_info.tsv", "w") as out:
        out.write("gene_name\tChr\tStart\tEnd\tStrand\tLength\n")
        for g in genes:
            out.write("\t".join([g, *gene_info[g]]) + "\n")
    print(f"{a.batch}: {len(samples)} samples, {len(genes)} genes")


if __name__ == "__main__":
    main()
