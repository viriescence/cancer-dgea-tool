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
# sys is used to read optional project names from the command line
import sys

# Path provides cleaner and platform-independent file paths 
from pathlib import Path

# Matplotlib is used to create the volcano plots 
import matplotlib

# "Agg" allows matplotlib to create image files without requiring a 
# graphical density. This is useful for running the script outside a notebook.
matplotlib.use("Agg")  # draw to files only, so it works without a display
import matplotlib.pyplot as plt

# NumPy is used for numerical operations such as checking for
# finite values and calculating -log10 adjusted p-values.
import numpy as np

# pandas is used to read, manipulate, merge, and save tabular data
import pandas as pd

# for welch's independent sample t-test
from scipy.stats import ttest_ind

# Used to perform Benjamimi-Hochberg multiple-testing correction
from statsmodels.stats.multitest import multipletests

# -------------------
# File locations
# -------------------

# Raw count data used as input for DESeq2
SAMPLE_DIR = Path("data") / "sample"

# the analysis ready files from the R pre processing
PROCESSED_DIR = Path("data") / "processed"

# Ranked differential-expression tables are written and saved here
TABLES_DIR = Path("results") / "tables"

# Volcano plots are written and saved here
FIGURES_DIR = Path("results") / "figures"

# -------------------
# File locations
# -------------------

# A gene must have an adjusted p-value below 0.05 to be
# considered statistically significant.
FDR_THRESHOLD = 0.05

# |log2FC| >= 1 corresponds to at least a two-fold difference
# in expression between tumor and normal samples.
LOG2FC_THRESHOLD = 1.0

# -------------------
# Load input data
# -------------------
def load_inputs(project, data_dir=PROCESSED_DIR):
    """
    Load the processed expression data, sample group labels, 
    and gene annotations for one TCGA project.

    Parameters 
    -----------
    project : str
        TCGA project name, such as "TCGA-BRCA" or "TCGA-LUAD"

    data_dir : Path
        This is a directory containing the processed input files

    Returns
    ----------
    expr : pandas.DataFrame
        Gene expression matrix with genes as rows and samples as columns

    group : panda.Series
        Tumor/normal label for each sample, ordered to match the columns of the expression matrix.

    genes : pandas.DataFrame
        Gene annotation label containing gene ID's and names
    """

    # Load the log2-CPM expression matrix.
    # Rows = genes 
    # Columns = patient/sample ID's
    expr = pd.read_csv(data_dir / f"{project}_logcpm.csv.gz", index_col=0)

    # Load the tumor/normal classification for each sample
    # The sample ID is used as a DataFrame index.
    # Columns = patient/sample IDs
    groups = pd.read_csv(data_dir / f"{project}_groups.csv", index_col=0)

     # Load gene annotations, including readable gene names.
    genes = pd.read_csv(data_dir / f"{project}_genes.csv")

    # Make sure sample IDs use the same string data type in 
    # both the expression and group tables
    expr.columns = expr.columns.astype(str)
    groups.index = groups.index.astype(str)

    # Find any expression samples that do not have a corresponding tumor/label
    missing = set(expr.columns) - set(groups.index)

    # Stop the analysis if sample labels are missing
    # Without a group label, we cannot determine whether a 
    # Sample belongs to a the tumor or normal group
    if missing:
        raise ValueError(f"Missing group labels for samples: {sorted(missing)[:10]}")

    # Reorder the group labels so they exactly match the column order of the expression matrix
    return expr, groups.loc[expr.columns, "group"], genes


# ----------------------
# Differential expression analysis
# ----------------------
def analyze_project(project, data_dir=PROCESSED_DIR):
    """
    Perform differential gene expression analysis for one
    TCGA cancer project.

    Tumor and normal expression values are compared against 
    Welch's t-test. P-values are corrected using the Benjamini-Hochberg
    false discovery rate procedure.

    Genes are classified as: 
        up   = adjusted p < 0.05 and log2FC >= 1
        down = adjusted p > 0.05 and log2FC <= -1
        ns   = does not meet both criteria
    """

    # -------------
    # Load the data
    # -------------
    expr, group, genes = load_inputs(project, data_dir)

    # -------------
    # Separate the sample groups 
    # -------------

    # Select expression columns belonging to normal samples
    normal = expr.loc[:, (group == "normal").to_numpy()]

    # Select expression columns belonging to tumor samples
    tumor = expr.loc[:, (group == "tumor").to_numpy()]

    # Display basic information about the dataset
    print(f"{project}: {normal.shape[1]} normal, {tumor.shape[1]} tumor, "
          f"{expr.shape[0]} genes")

    # -------------
    # Calculate log2 fold change 
    # -------------

    # Calculate the average expression for every gene in normal and tumor samples
    mean_normal = normal.mean(axis=1)
    mean_tumor = tumor.mean(axis=1)

    # The expression data is already in log2 scale, 
    # so mean tumor and mean normal can be directly subtracted
    # Positive values indicate high expression, negative indicates lower
    log2fc = mean_tumor - mean_normal

    # --------------
    # Welch's t-test
    # --------------

    # Perform an independent t-test per gene
    # equal_var = False tells scipy not to assume that the tumor and normal variants have equal variances
    # axis = 1 means each row(gene) receives its own test
    t_statistic, p_value = ttest_ind(
        tumor.to_numpy(dtype=float),
        normal.to_numpy(dtype=float),
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )

    # If a NaN or infinite values are returned, replace the p-value with p = 1
    # They'll get treated as nonsignificant rather than breaking later calculations
    p_value = np.nan_to_num(p_value, nan=1.0, posinf=1.0, neginf=1.0)

    # -------------
    # Multiple testing correction
    # -------------

    # Thousands of genes are tested simultaneously.
    # Testing this many hypotheses increases the probability 
    # of obtaining apparently significant results
    # Benjamin Hochberg (BH) correction controls the false
    # discovery rates (FDR)
    padj = multipletests(p_value, method="fdr_bh")[1]

    # -------------
    # Build results table
    # -------------

    # Create one more row per gene containing its fold change and
    # statistical-test results
    results = pd.DataFrame({
        "gene_id": expr.index.astype(str),
        "log2FC": log2fc.to_numpy(),
        "p_value": p_value,
        "padj": padj,
    })

    # Convert gene IDs to strings so that they can be matched
    # reliably between the results and annotation tables
    genes["gene_id"] = genes["gene_id"].astype(str)

    # Add human readable gene names to the statistical results
    results = results.merge(
        genes[["gene_id", "gene_name"]].drop_duplicates("gene_id"),
        on="gene_id",
        how="left",
    )

    # ----------
    # Classify differential expression
    # ----------

    # Determine which genes pass the FDR threshold
    significant = results["padj"] < FDR_THRESHOLD

    # Start by classifying every gene as nonsignificant.
    results["status"] = "ns"

    # Upregulated if statistically significant AND log2FC >= +1.
    results.loc[significant & (results["log2FC"] >= LOG2FC_THRESHOLD), "status"] = "up"

    # Downregulated if statistically significant and log2FC <= -1
    results.loc[significant & (results["log2FC"] <= -LOG2FC_THRESHOLD), "status"] = "down"

    # -------------
    # Rank genes
    # -------------

    # Temporarily calculate absolute fold change so that
    # genes with equally significant adjusted p-values can be ranked by magnitude of their expression difference
    results = (
        results.assign(abs_log2FC=results["log2FC"].abs())
        .sort_values(["padj", "abs_log2FC"], ascending=[True, False])
        .drop(columns="abs_log2FC")
        .reset_index(drop=True)
    )

    # Keep the final columns in a convenient order
    results = results[["gene_id", "gene_name", "log2FC", "p_value", "padj", "status"]]

    # ------------
    # Save the ranked gene table
    # ------------

    # Create the output directory if it does not already exist
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # Construct the output filename
    table_path = TABLES_DIR / f"{project}_ranked_genes.csv"

    # Save all ranked genes to CSV
    results.to_csv(table_path, index=False)

    # Quick summary for the console: counts per status and the top 10 genes.
    print(results["status"].value_counts().to_string())

    # Print a quick summary to the console 
    print(results.head(10).to_string(index=False))

    # Display the 10 highest ranked genes 
    print(f"Saved {table_path}")

    # Create the volcano plot 
    figure_path = save_volcano(results, project)
    print(f"Saved {figure_path}")

    # Return the full results so it can be used in python or a notebook
    return results


def save_volcano(results, project):
    """
    Create and save a volcano plot
    X-axis:
        log2 fold change between tumor and normal samples
    Y-axis:
        -log10 adjusted p-value

    Colors: 
        Gray = not significant 
        Red = upregulated
        Blue = downregulated
    """
    # Remove genes with unusable fold change or adjusted
    # p values before plotting
    plot_data = results[
        np.isfinite(results["log2FC"]) &
        np.isfinite(results["padj"])
    ].copy()

    # The floor keeps -log10 finite when a p-value underflows to zero.
    # Smaller adjusted p-values therefore appear higher on volcano plots
    # clip() prevents log(10), which would produce infinity
    plot_data["neg_log10_padj"] = -np.log10(
        plot_data["padj"].clip(lower=np.finfo(float).tiny)
    )

    # Give a color to each differential expression category
    colors = {
        "ns": "lightgray",
        "up": "firebrick",
        "down": "royalblue"
    }

    # Create the plotting area
    fig, ax = plt.subplots(figsize=(8, 6))

    # One scatter layer per status so each gets its own colour and legend entry.
    for status, color in colors.items():
        subset = plot_data[plot_data["status"] == status]

        ax.scatter(
            subset["log2FC"],
            subset["neg_log10_padj"],
            color=color,
            s=10,
            alpha=0.6,
            label=status
        )

    

    # -----------
    # Draw Significance thresholds
    # -----------
    
    # Dashed lines mark the two thresholds betwee log2FC = -1 and +1
    #
    # These correspond to the two fold expression change used to classify up/downregulated genes
    for x in (-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD):
        ax.axvline(
            x,
            color="black",
            linestyle="--",
            linewidth=1
        )
    # Horizontal line marks adjusted p = 0.05
    ax.axhline(
        -np.log10(FDR_THRESHOLD),
        color="black",
        linestyle="--",
        linewidth=1
    )

    # -------------
    # Label top 10 genes
    # -------------

    # Label the 10 most significant up/down genes (the table is sorted by padj).
    # Keep only genes classified as significantly up or downregulated
    # Because the results table was previously sorted by adjusted p value, head(10) 
    # should select the 10 most significantly genes
    top_genes = plot_data[
        plot_data["status"].isin(["up", "down"])
    ].head(10)

    # Add the gene name besides each selected point 
    for _, gene in top_genes.iterrows():
        # Use the human readable gene symbol when available
        label = gene["gene_name"]

        # Fall back to the ensemble gene ID if no gene name is available
        if pd.isna(label):
            label = gene["gene_id"]

        ax.annotate(
            label,
            # location of the gene's point
            xy=(
                gene["log2FC"],
                gene["neg_log10_padj"]
            ),
            # Move text away from the point
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
            # Draw a small line connecting the label to
            # its corresponding gene.
            arrowprops={
                "arrowstyle": "-",
                "linewidth": 0.5
            }
        )

    # ----------
    # Plot formatting
    # ----------
    ax.set_title(
        f"{project} - Python (Welch t-test, BH)"
    )
    ax.set_xlabel(
        "log2 fold change (tumor vs. normal)"
    )
    ax.set_ylabel(
        "-log10 adjusted p-value"
    )

    # display the ns/up/down legend
    ax.legend()

    # automatically adjust spacing so labels and axes fit 
    fig.tight_layout()

    # -----------------------------
    # Save volcano plot
    # -----------------------------
    
    # Create the figure directory if necessary
    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # Construct the output filename
    path = FIGURES_DIR / f"{project}_volcano_plot.png"

    # Save the figure as a PNG
    fig.savefig(
        path,
        dpi=150,
        bbox_inches="tight"
    )

    # close the matplotlib figure
    plt.close(fig)

    # if ran on colab or jupyter, display the notebook
    show_in_notebook(path)

    
    return path


def compare_with_deseq2(welch_results, project):
    """
    Compare Welch t-test with DESeq2
    Compares: 
        1. Pearson correlation of log2 fold change
        2. Spearmen correlation of log2 fold changes
        3. Agreement and fold-change direction 
        4. Agreement in up/down/ns classification
        5. Scatterplot of DESeq2 vs Welch log2 fold changes
    """
    # -------------
    # Load DESeq2 results
    # ------------

    # Change this path if your DESeq2 results are stored elsewhere
    deseq2_path = (
        TABLES_DIR /
        f"{project}_deseq2.csv"
    )

    deseq2 = pd.read_csv(deseq2_path)

    print("DESeq2 columns:")
    print(deseq2.columns.tolist())
    print(deseq2.head())

    # Make sure gene IDs have the same data type.
    deseq2["gene_id"] = deseq2["gene_id"].astype(str)

    welch_results = welch_results.copy()
    welch_results["gene_id"] = (
        welch_results["gene_id"].astype(str)
    )


    # ----------------
    # Prepare DESeq2 columns
    # ----------------

    # DESeq2 normally calls these columns:
    #
    # log2FoldChange
    # pvalue
    # padj
    #
    # Rename them so it is obvious which method they came from.
    deseq2 = deseq2.rename(
        columns={
            "log2FoldChange": "log2FC_deseq2",
            "pvalue": "p_value_deseq2",
            "padj": "padj_deseq2"
        }
    )

    # Rename the corresponding Welch columns.
    welch = welch_results.rename(
        columns={
            "log2FC": "log2FC_welch",
            "p_value": "p_value_welch",
            "padj": "padj_welch",
            "status": "status_welch"
        }
    )


    # -----------------------------
    # Merge results by gene
    # -----------------------------

    # Only genes found in both analyses are compared.
    comparison = welch.merge(
        deseq2[
            [
                "gene_id",
                "log2FC_deseq2",
                "p_value_deseq2",
                "padj_deseq2"
            ]
        ],
        on="gene_id",
        how="inner"
    )

    # Remove genes without usable fold-change estimates.
    comparison = comparison[
        np.isfinite(comparison["log2FC_welch"]) &
        np.isfinite(comparison["log2FC_deseq2"])
    ].copy()

    print(
        f"\n{project}: "
        f"{len(comparison)} genes compared"
    )


    # -----------------------------
    # Pearson correlation
    # -----------------------------

    # Pearson measures how closely the numerical fold-change
    # estimates follow a linear relationship.
    pearson_r, pearson_p = pearsonr(
        comparison["log2FC_deseq2"],
        comparison["log2FC_welch"]
    )

    print(
        f"Pearson correlation: "
        f"r = {pearson_r:.4f}"
    )


    # -----------------------------
    # Spearman correlation
    # -----------------------------

    # Spearman measures whether the two methods rank genes
    # similarly, without requiring a perfectly linear relationship.
    spearman_rho, spearman_p = spearmanr(
        comparison["log2FC_deseq2"],
        comparison["log2FC_welch"]
    )

    print(
        f"Spearman correlation: "
        f"rho = {spearman_rho:.4f}"
    )


    # -----------------------------
    # Direction agreement
    # -----------------------------

    # Determine whether both methods give the gene the same
    # direction of fold change.
    #
    # Positive = higher in tumor
    # Negative = lower in tumor
    comparison["direction_agrees"] = (
        np.sign(comparison["log2FC_welch"]) ==
        np.sign(comparison["log2FC_deseq2"])
    )

    direction_agreement = (
        comparison["direction_agrees"].mean() * 100
    )

    print(
        f"Direction agreement: "
        f"{direction_agreement:.2f}%"
    )


    # -----------------------------
    # Classify DESeq2 results
    # -----------------------------

    # Start every DESeq2 gene as nonsignificant.
    comparison["status_deseq2"] = "ns"

    deseq2_significant = (
        comparison["padj_deseq2"] < FDR_THRESHOLD
    )

    # Significant and at least two-fold higher in tumor.
    comparison.loc[
        deseq2_significant &
        (
            comparison["log2FC_deseq2"]
            >= LOG2FC_THRESHOLD
        ),
        "status_deseq2"
    ] = "up"

    # Significant and at least two-fold lower in tumor.
    comparison.loc[
        deseq2_significant &
        (
            comparison["log2FC_deseq2"]
            <= -LOG2FC_THRESHOLD
        ),
        "status_deseq2"
    ] = "down"


    # -----------------------------
    # Compare classifications
    # -----------------------------

    # Create a table showing how often Welch and DESeq2
    # agree on up/down/not-significant classifications.
    agreement_table = pd.crosstab(
        comparison["status_welch"],
        comparison["status_deseq2"],
        rownames=["Welch"],
        colnames=["DESeq2"]
    )

    print("\nClassification agreement:")
    print(agreement_table)


    # Calculate the overall percentage of identical classifications.
    comparison["classification_agrees"] = (
        comparison["status_welch"] ==
        comparison["status_deseq2"]
    )

    classification_agreement = (
        comparison["classification_agrees"].mean()
        * 100
    )

    print(
        f"\nOverall classification agreement: "
        f"{classification_agreement:.2f}%"
    )


    # -----------------------------
    # Save comparison table
    # -----------------------------

    comparison_path = (
        TABLES_DIR /
        f"{project}_Welch_vs_DESeq2.csv"
    )

    comparison.to_csv(
        comparison_path,
        index=False
    )

    print(
        f"Saved comparison table to: "
        f"{comparison_path}"
    )


    # -----------------------------
    # Fold-change comparison plot
    # -----------------------------

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    ax.scatter(
        comparison["log2FC_deseq2"],
        comparison["log2FC_welch"],
        s=10,
        alpha=0.4
    )

    # Determine the common plotting range.
    low = min(
        comparison["log2FC_deseq2"].min(),
        comparison["log2FC_welch"].min()
    )

    high = max(
        comparison["log2FC_deseq2"].max(),
        comparison["log2FC_welch"].max()
    )

    # Perfect agreement would fall on this diagonal line.
    ax.plot(
        [low, high],
        [low, high],
        linestyle="--",
        linewidth=1
    )

    # Add zero lines to make disagreements in direction
    # easier to see.
    ax.axhline(
        0,
        linestyle=":",
        linewidth=1
    )

    ax.axvline(
        0,
        linestyle=":",
        linewidth=1
    )

    ax.set_xlabel(
        "DESeq2 log2 fold change"
    )

    ax.set_ylabel(
        "Welch log2 fold change"
    )

    ax.set_title(
        f"{project}: DESeq2 vs Welch\n"
        f"Pearson r = {pearson_r:.3f}, "
        f"Spearman ρ = {spearman_rho:.3f}"
    )

    fig.tight_layout()


    # -----------------------------
    # Save comparison plot
    # -----------------------------

    FIGURES_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    plot_path = (
        FIGURES_DIR /
        f"{project}_Welch_vs_DESeq2.png"
    )

    fig.savefig(
        plot_path,
        dpi=150,
        bbox_inches="tight"
    )

    plt.close(fig)

    show_in_notebook(plot_path)

    print(
        f"Saved comparison plot to: "
        f"{plot_path}"
    )

    return comparison, agreement_table
def show_in_notebook(path):
    """
    Display the saved plot when running in Colab or Jupyter.
    When ran from terminal, this won't do anything
    """
    # IPython is usually available with jupyter or colab, but can't run from terminal
    try:
        from IPython import get_ipython
        from IPython.display import Image, display
    except ImportError:
        return
    # only displays the image when an IPython environment is actually active
    if get_ipython() is not None:
        display(Image(filename=str(path)))

# ------------
# Run analysis 
# ------------
if __name__ == "__main__":
    # Notebook kernels (Colab, Jupyter) add their own flags; ignore those.
    projects = [a for a in sys.argv[1:]
                if not a.startswith("-") and not a.endswith(".json")]
    # If no projects were explicitly provided, analyze breast cancer (BRCA) and Lung cancer (LUAD)
    projects = projects or ["TCGA-BRCA", "TCGA-LUAD"]
    # Run the analysis independently for each project
    for project in projects:
        welch_results = analyze_project(project)

        # Compare welch_results against DESeq2
        comparison, agreemnt_table = compare_with_deseq2(welch_results, project)
