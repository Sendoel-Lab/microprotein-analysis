"""Loader for the per-cluster glmGamPoi differential-expression tables of the
in vivo single-cell CRISPR screen (all tested gene x perturbation pairs).

Expected files: ${MP_EXTERNAL}/perturbseq/glm/p60_ko_vs_ctrl_cluster_<i>_glm.csv.gz
with columns name, pval, adj_pval, lfc, perturb.
"""
import os
from pathlib import Path

import pandas as pd

GLM_DIR = Path(os.environ["MP_EXTERNAL"]) / "perturbseq" / "glm"
GLM_PATTERN = "p60_ko_vs_ctrl_cluster_{}_glm.csv.gz"
CLUSTERS = range(11)


def load_unfiltered_glm(clusters=CLUSTERS, verbose=True):
    """Return {'Cluster i': DataFrame(name, pval, adj_pval, lfc, perturb, cluster)}."""
    out = {}
    for i in clusters:
        label = f"Cluster {i}"
        df = pd.read_csv(
            GLM_DIR / GLM_PATTERN.format(i),
            usecols=["name", "pval", "adj_pval", "lfc", "perturb"],
            dtype={"name": "category", "perturb": "category", "lfc": "float32"},
        )
        df["cluster"] = label
        out[label] = df
        if verbose:
            print(f"  {label}: {len(df):,} rows, {df['perturb'].nunique()} perturbations, "
                  f"{df['name'].nunique()} genes")
    return out
