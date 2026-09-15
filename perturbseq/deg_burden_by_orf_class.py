#!/usr/bin/env python3
"""Number and fraction of differentially expressed genes per perturbation,
compared between ORF classes (uORF, novel ORF, main ORF, 5'UTR), aggregated over
all clusters and per cluster.

Usage: perturbseq/deg_burden_by_orf_class.py
Inputs:  ${MP_EXTERNAL}/perturbseq/sgRNAs_ORFtype.xlsx (columns sgrna, ORF_type)
         glmGamPoi tables via glm_loader
Outputs: ${MP_RESULTS}/perturbseq/deg_burden/
         perturbation_deg_summary_<view>.csv, deg_burden_stats_<view>.txt,
         deg_burden_by_orf_class_<view>.{svg,png}
"""
import os
import re
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from glm_loader import load_unfiltered_glm

warnings.filterwarnings("ignore", category=FutureWarning)
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 10, "axes.linewidth": 1.0, "xtick.major.width": 0.8,
                     "ytick.major.width": 0.8, "svg.fonttype": "none"})

SGRNA_ANNOT = Path(os.environ["MP_EXTERNAL"]) / "perturbseq" / "sgRNAs_ORFtype.xlsx"
RESULTS = Path(os.environ["MP_RESULTS"]) / "perturbseq" / "deg_burden"
RESULTS.mkdir(parents=True, exist_ok=True)

PADJ, LFC = 0.05, 0.25
ORF_ORDER = ["upstream ORF", "new ORF", "main ORF", "5UTR"]
ORF_LABELS = ["uORF", "novel ORF", "main ORF", "5'UTR"]
ORF_COLORS = {"upstream ORF": "#D64550", "new ORF": "#3B8EA5", "main ORF": "#6B9F78", "5UTR": "#E8985E"}

## ---- perturbation -> ORF class ----
sg = pd.read_excel(SGRNA_ANNOT)
sg["perturb_label"] = sg["sgrna"].str.replace(r"_sgRNA\.\d+$", "", regex=True)
sg["orf_type_simple"] = sg["ORF_type"].replace({"upstream ORF1": "upstream ORF", "upstream ORF2": "upstream ORF"})
perturb_to_type = sg.groupby("perturb_label")["orf_type_simple"].first().to_dict()
perturb_to_gene = {p: re.sub(r"_(mORF|uORF\d*|5UTR|\d+aa)$", "", p) for p in perturb_to_type}

## ---- DE tables ----
cluster_all = load_unfiltered_glm()
cluster_degs = {k: df[(df["adj_pval"] < PADJ) & (df["lfc"].abs() > LFC)].copy() for k, df in cluster_all.items()}


def build_summary(all_df, degs_df, view):
    n_tested = all_df.groupby("perturb")["name"].nunique().rename("n_tested")
    agg = degs_df.groupby("perturb").agg(
        n_deg=("name", "nunique"),
        n_deg_entries=("name", "count"),
        mean_abs_lfc=("lfc", lambda x: np.mean(np.abs(x[x.abs() < 100]))),
        clusters_affected=("cluster", "nunique"))
    s = n_tested.to_frame().join(agg, how="left")
    s[["n_deg", "n_deg_entries", "clusters_affected"]] = s[["n_deg", "n_deg_entries", "clusters_affected"]].fillna(0).astype(int)
    s["mean_abs_lfc"] = s["mean_abs_lfc"].fillna(0.0)
    s["deg_fraction"] = s["n_deg"] / s["n_tested"]
    s = s.reset_index()
    s["orf_type"] = s["perturb"].map(perturb_to_type)
    s["host_gene"] = s["perturb"].map(perturb_to_gene)
    s = s[s["orf_type"].notna() & (s["orf_type"] != "non targeting control")]
    s.to_csv(RESULTS / f"perturbation_deg_summary_{view}.csv", index=False)
    return s


def run_stats(s, view):
    uorf = s[s["orf_type"] == "upstream ORF"]
    new = s[s["orf_type"] == "new ORF"]
    if len(uorf) < 2 or len(new) < 2:
        return None
    out = [f"DEG burden ({view}): upstream ORF vs new ORF"]
    for label, g in [("upstream ORF", uorf), ("new ORF", new)]:
        out.append(f"{label} (n={len(g)}): n_deg median={g['n_deg'].median():.0f}, "
                   f"deg_fraction median={g['deg_fraction'].median():.4f}, "
                   f"zero-DEG perturbations={(g['n_deg'] == 0).sum()}")
    res = {}
    for metric in ("n_deg", "deg_fraction", "n_tested"):
        u, p = stats.mannwhitneyu(uorf[metric], new[metric], alternative="two-sided")
        out.append(f"Mann-Whitney U on {metric}: U={u:.0f}, p={p:.2e}")
        res[metric] = p

    both = s[s["orf_type"].isin(["upstream ORF", "new ORF"])].copy()
    both["is_uorf"] = both["orf_type"] == "upstream ORF"
    for metric in ("n_deg", "deg_fraction"):
        q75 = s[metric].quantile(0.75)
        ct = pd.crosstab(both["is_uorf"], both[metric] >= q75)
        if ct.shape == (2, 2):
            odds, p = stats.fisher_exact(ct)
            out.append(f"Fisher's exact, uORF in top quartile of {metric} (>= {q75:.4g}): OR={odds:.2f}, p={p:.4f}")

    groups = [g["deg_fraction"].values for _, g in s.groupby("orf_type") if len(g) >= 2]
    if len(groups) >= 2:
        h, p = stats.kruskal(*groups)
        out.append(f"Kruskal-Wallis on deg_fraction across ORF classes: H={h:.2f}, p={p:.2e}")

    text = "\n".join(out)
    (RESULTS / f"deg_burden_stats_{view}.txt").write_text(text + "\n")
    print(text)
    return res


def boxstrip(ax, d, metric, pval):
    data = [d[d["orf_type"] == o][metric].values for o in ORF_ORDER]
    idx = [i for i, v in enumerate(data) if len(v) > 0]
    bp = ax.boxplot([data[i] for i in idx], positions=idx, widths=0.5, patch_artist=True,
                    showfliers=False, zorder=2, medianprops=dict(color="white", linewidth=2),
                    whiskerprops=dict(color="#555555", linewidth=1), capprops=dict(color="#555555", linewidth=1))
    for patch, i in zip(bp["boxes"], idx):
        patch.set(facecolor=ORF_COLORS[ORF_ORDER[i]], alpha=0.7, edgecolor="white", linewidth=1)
    rng = np.random.default_rng(42)
    for i in idx:
        ax.scatter(i + rng.normal(0, 0.08, len(data[i])), data[i], color=ORF_COLORS[ORF_ORDER[i]],
                   s=18, alpha=0.5, edgecolors="white", linewidths=0.3, zorder=3)
        ax.text(i, -0.08, f"n={len(data[i])}", transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=9, color="#666666")
    ax.set_xticks(range(len(ORF_ORDER)))
    ax.set_xticklabels(ORF_LABELS, fontsize=11)
    if 0 in idx and 1 in idx:
        y = max(data[0].max(), data[1].max()) * 1.3 or 1.0
        ax.plot([0, 0, 1, 1], [y * 0.85, y, y, y * 0.85], lw=1.2, color="black", clip_on=False)
        ax.text(0.5, y * 1.08, f"p = {pval:.1e}" if pval < 0.001 else f"p = {pval:.3f}",
                ha="center", va="bottom", fontsize=9, fontstyle="italic")


def make_figure(s, view, title, res):
    d = s[s["orf_type"].isin(ORF_ORDER)].copy()
    d["orf_type"] = pd.Categorical(d["orf_type"], categories=ORF_ORDER, ordered=True)
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35, height_ratios=[1.2, 1])
    fig.suptitle(title, fontsize=13, color="#555555", y=0.98)

    ax = fig.add_subplot(gs[0, 0])
    boxstrip(ax, d, "n_deg", res["n_deg"])
    ax.set_yscale("symlog", linthresh=10)
    ax.set_ylabel(f"Number of DEGs\n(|log2FC| > {LFC}, adj. p < {PADJ})", fontsize=11)
    ax.set_title("A", fontsize=14, fontweight="bold", loc="left", pad=8)
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[0, 1])
    boxstrip(ax, d, "deg_fraction", res["deg_fraction"])
    ax.set_ylabel("DEG fraction\n(n_deg / n_tested)", fontsize=11)
    ax.set_title("B", fontsize=14, fontweight="bold", loc="left", pad=8)
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 0])
    for orf, label in zip(ORF_ORDER[:2], ORF_LABELS[:2]):
        vals = np.sort(d[d["orf_type"] == orf]["deg_fraction"].values)
        if len(vals):
            ax.step(vals, 1 - np.arange(1, len(vals) + 1) / len(vals), where="post",
                    color=ORF_COLORS[orf], linewidth=2, label=label)
    ax.set_xlabel("DEG fraction", fontsize=11)
    ax.set_ylabel("Fraction of perturbations\nwith >= x", fontsize=11)
    ax.axhline(0.25, ls=":", color="gray", alpha=0.5, lw=0.8)
    ax.legend(frameon=False, fontsize=10, loc="upper right")
    ax.set_title("C", fontsize=14, fontweight="bold", loc="left", pad=8)
    sns.despine(ax=ax)

    ax = fig.add_subplot(gs[1, 1])
    ranked = d.sort_values("n_deg", ascending=False).reset_index(drop=True)
    ax.bar(range(len(ranked)), ranked["n_deg"].values, color=[ORF_COLORS[t] for t in ranked["orf_type"]],
           width=1.0, edgecolor="none", alpha=0.8)
    ax.set_xlabel("Perturbations (ranked)", fontsize=11)
    ax.set_ylabel("Number of DEGs", fontsize=11)
    ax.legend(handles=[Patch(facecolor=ORF_COLORS[o], label=l, alpha=0.8) for o, l in zip(ORF_ORDER, ORF_LABELS)],
              frameon=False, fontsize=9, ncol=2, loc="upper right")
    ax.set_title("D", fontsize=14, fontweight="bold", loc="left", pad=8)
    sns.despine(ax=ax)

    stem = RESULTS / f"deg_burden_by_orf_class_{view}"
    plt.savefig(f"{stem}.svg", bbox_inches="tight", facecolor="white")
    plt.savefig(f"{stem}.png", dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()


views = {"all_clusters": (pd.concat(cluster_all.values(), ignore_index=True),
                          pd.concat(cluster_degs.values(), ignore_index=True), "All clusters")}
for label in cluster_all:
    views[label.replace("Cluster ", "cluster")] = (cluster_all[label], cluster_degs[label], label)

for view, (all_df, degs_df, title) in views.items():
    if len(degs_df) < 50:
        continue
    print(f"\n== {title} ==")
    summary = build_summary(all_df, degs_df, view)
    res = run_stats(summary, view)
    if res:
        make_figure(summary, view, title, res)
