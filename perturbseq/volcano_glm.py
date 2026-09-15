#!/usr/bin/env python3
"""Volcano plots of glmGamPoi differential expression for selected perturbations
within one cell cluster.

Usage: perturbseq/volcano_glm.py [--cluster 0] [--perturbations Emg1 Gm13935] [--lfc 0.5]
Output: ${MP_RESULTS}/perturbseq/volcano/volcano_<perturbation>_lfc<threshold>.{pdf,png}
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                     "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"]})
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from adjustText import adjust_text

from glm_loader import load_unfiltered_glm

PADJ = 0.05
PADJ_FLOOR = 1e-300
N_LABELS = 24
UP, DOWN, NS = "firebrick", "royalblue", "lightgrey"
PLAIN = FuncFormatter(lambda v, _: f"{v:g}")


def volcano(df, gene, lfc_thr, cluster, out_dir):
    d = df[df["perturb"] == gene].copy()
    d["nlp"] = -np.log10(d["adj_pval"].clip(lower=PADJ_FLOOR))
    sig = d["adj_pval"] < PADJ
    up, down = sig & (d["lfc"] >= lfc_thr), sig & (d["lfc"] <= -lfc_thr)
    ns = ~(up | down)

    fig, ax = plt.subplots(figsize=(5.0, 5.4))
    ax.scatter(d.loc[ns, "lfc"], d.loc[ns, "nlp"], s=6, c=NS, linewidths=0, rasterized=True)
    for m, col in [(down, DOWN), (up, UP)]:
        ax.scatter(d.loc[m, "lfc"], d.loc[m, "nlp"], s=7, c=col, linewidths=0)
    for x in (lfc_thr, -lfc_thr):
        ax.axvline(x, color="0.4", linestyle=":", linewidth=0.9)
    ax.axhline(-np.log10(PADJ), color="0.4", linestyle=":", linewidth=0.9)

    lab = d[up | down].nlargest(N_LABELS, "nlp")
    texts = [ax.text(r.lfc, r.nlp, r["name"], fontsize=7, fontstyle="italic") for _, r in lab.iterrows()]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="0.6", lw=0.4), expand=(1.2, 1.4))

    xmax = np.ceil(d["lfc"].abs().max() * 10) / 10
    ymax = d["nlp"].max() * 1.12
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(-2, ymax)
    ax_x = (d["lfc"] + xmax) / (2 * xmax)
    in_top = (d["nlp"] + 2) / (ymax + 2) > 0.82
    left_free = int((in_top & (ax_x < 0.42)).sum()) <= int((in_top & (ax_x > 0.58)).sum())
    note_x, note_ha = (0.02, "left") if left_free else (0.98, "right")
    ax.text(note_x, 0.98, f"glmGamPoi\nadj p < {PADJ}, |log2FC| >= {lfc_thr:g}",
            transform=ax.transAxes, ha=note_ha, va="top", fontsize=6.5, color="0.3", linespacing=1.4)

    ax.set_xlabel(f"log2 Fold Change Transcriptome ({gene} vs. control)")
    ax.set_ylabel("-log10(p adj)")
    ax.set_title(f"{gene} perturbation, cluster {cluster}")
    ax.xaxis.set_major_formatter(PLAIN)
    ax.yaxis.set_major_formatter(PLAIN)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()

    stem = out_dir / f"volcano_{gene}_lfc{lfc_thr:g}"
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png", dpi=300)
    plt.close(fig)
    print(f"{gene} |lfc| >= {lfc_thr:g}: up = {int(up.sum())}, down = {int(down.sum())}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cluster", type=int, default=0)
    ap.add_argument("--perturbations", nargs="+", default=["Emg1", "Gm13935"])
    ap.add_argument("--lfc", type=float, nargs="+", default=[0.5])
    a = ap.parse_args()
    out_dir = Path(os.environ["MP_RESULTS"]) / "perturbseq" / "volcano"
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_unfiltered_glm([a.cluster], verbose=False)[f"Cluster {a.cluster}"]
    df = df[df["lfc"].abs() < 50]
    for gene in a.perturbations:
        for thr in a.lfc:
            volcano(df, gene, thr, a.cluster, out_dir)


if __name__ == "__main__":
    main()
