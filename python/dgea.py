"""Differential gene expression of tumor vs. normal samples (Python).

Run from the project root, after `Rscript run_all.R` has written the inputs:

    python3 python/dgea.py                  # TCGA-BRCA and TCGA-LUAD
    python3 python/dgea.py TCGA-COAD ...    # other projects

Reads, from data/processed/:
    <project>_logcpm.csv.gz   genes x samples, log2 CPM
    <project>_groups.csv      sample_id, group ("tumor" / "normal")
    <project>_genes.csv       gene_id, gene_name, gene_type

For each gene: log2 fold change (mean tumor - mean normal on log2 CPM),
Welch's t-test, Benjamini-Hochberg adjusted p-value. A gene is "up" or "down"
at adjusted p < 0.05 and |log2 fold change| >= 1, otherwise "ns".

Writes:
    results/tables/<project>_ranked_genes.csv
    results/figures/<project>_volcano_plot.png
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw to files only, so it works without a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

PROCESSED_DIR = Path("data") / "processed"
TABLES_DIR = Path("results") / "tables"
FIGURES_DIR = Path("results") / "figures"

FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0


def load_inputs(project, data_dir=PROCESSED_DIR):
    """Read the three files R exported and line up the group labels."""
    expr = pd.read_csv(data_dir / f"{project}_logcpm.csv.gz", index_col=0)
    groups = pd.read_csv(data_dir / f"{project}_groups.csv", index_col=0)
    genes = pd.read_csv(data_dir / f"{project}_genes.csv")

    expr.columns = expr.columns.astype(str)
    groups.index = groups.index.astype(str)

    # Every sample in the expression table must have a tumor/normal label.
    missing = set(expr.columns) - set(groups.index)
    if missing:
        raise ValueError(f"Missing group labels for samples: {sorted(missing)[:10]}")

    # Put the labels in the same order as the expression columns.
    return expr, groups.loc[expr.columns, "group"], genes


def analyze_project(project, data_dir=PROCESSED_DIR):
    """Run the full analysis for one project and save its table and plot."""
    expr, group, genes = load_inputs(project, data_dir)

    # Split the columns into the two groups being compared.
    normal = expr.loc[:, (group == "normal").to_numpy()]
    tumor = expr.loc[:, (group == "tumor").to_numpy()]
    print(f"{project}: {normal.shape[1]} normal, {tumor.shape[1]} tumor, "
          f"{expr.shape[0]} genes")

    # The data is already log2, so a difference of means is a log2 fold change.
    mean_normal = normal.mean(axis=1)
    mean_tumor = tumor.mean(axis=1)
    log2fc = mean_tumor - mean_normal

    # One Welch t-test per gene (equal_var=False: no equal-variance assumption).
    t_statistic, p_value = ttest_ind(
        tumor.to_numpy(dtype=float),
        normal.to_numpy(dtype=float),
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)

    # Thousands of tests give false hits by chance; Benjamini-Hochberg controls
    # the false discovery rate.
    padj = multipletests(p_value, method="fdr_bh")[1]

    results = pd.DataFrame({
        "gene_id": expr.index.astype(str),
        "log2FC": log2fc.to_numpy(),
        "p_value": p_value,
        "padj": padj,
    })

    # Add readable gene names.
    genes["gene_id"] = genes["gene_id"].astype(str)
    results = results.merge(
        genes[["gene_id", "gene_name"]].drop_duplicates("gene_id"),
        on="gene_id",
        how="left",
    )

    # A gene must be both significant and changed by at least two-fold.
    significant = results["padj"] < FDR_THRESHOLD
    results["status"] = "ns"
    results.loc[significant & (results["log2FC"] >= LOG2FC_THRESHOLD), "status"] = "up"
    results.loc[significant & (results["log2FC"] <= -LOG2FC_THRESHOLD), "status"] = "down"

    # Rank by adjusted p-value, breaking ties by size of the change.
    results = (
        results.assign(abs_log2FC=results["log2FC"].abs())
        .sort_values(["padj", "abs_log2FC"], ascending=[True, False])
        .drop(columns="abs_log2FC")
        .reset_index(drop=True)
    )
    results = results[["gene_id", "gene_name", "log2FC", "p_value", "padj", "status"]]

    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    table_path = TABLES_DIR / f"{project}_ranked_genes.csv"
    results.to_csv(table_path, index=False)

    # Quick summary for the console: counts per status and the top 10 genes.
    print(results["status"].value_counts().to_string())
    print(results.head(10).to_string(index=False))
    print(f"Saved {table_path}")

    figure_path = save_volcano(results, project)
    print(f"Saved {figure_path}")
    return results


def save_volcano(results, project):
    """Volcano plot: fold change on x, significance on y, coloured by status."""
    # Drop genes with no usable fold change or p-value.
    plot_data = results[np.isfinite(results["log2FC"]) & np.isfinite(results["padj"])].copy()
    # The floor keeps -log10 finite when a p-value underflows to zero.
    plot_data["neg_log10_padj"] = -np.log10(
        plot_data["padj"].clip(lower=np.finfo(float).tiny)
    )

    # One scatter layer per status so each gets its own colour and legend entry.
    colors = {"ns": "lightgray", "up": "firebrick", "down": "royalblue"}
    fig, ax = plt.subplots(figsize=(8, 6))
    for status, color in colors.items():
        subset = plot_data[plot_data["status"] == status]
        ax.scatter(subset["log2FC"], subset["neg_log10_padj"],
                   color=color, s=10, alpha=0.6, label=status)

    # Dashed lines mark the two thresholds.
    for x in (-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD):
        ax.axvline(x, color="black", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(FDR_THRESHOLD), color="black", linestyle="--", linewidth=1)

    ax.set_title(f"{project} - Python (Welch t-test, BH)")
    ax.set_xlabel("log2 fold change (tumor vs. normal)")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.legend()
    fig.tight_layout()

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    path = FIGURES_DIR / f"{project}_volcano_plot.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    show_in_notebook(path)
    return path


def show_in_notebook(path):
    """Display the saved plot when running in Colab or Jupyter; no-op otherwise."""
    try:
        from IPython import get_ipython
        from IPython.display import Image, display
    except ImportError:
        return
    if get_ipython() is not None:
        display(Image(filename=str(path)))


if __name__ == "__main__":
    # Notebook kernels (Colab, Jupyter) add their own flags; ignore those.
    projects = [a for a in sys.argv[1:]
                if not a.startswith("-") and not a.endswith(".json")]
    projects = projects or ["TCGA-BRCA", "TCGA-LUAD"]
    for project in projects:
        analyze_project(project)
