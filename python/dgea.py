"""dgea.py -- Analysis side (Daniel), Python version.

Reads the files written by export_for_python() in R/load_data.R:
    data/processed/<project>_logcpm.csv.gz   genes x samples, log2 CPM
    data/processed/<project>_groups.csv      sample_id, group ("tumor"/"normal")
    data/processed/<project>_genes.csv       gene_id, gene_name, gene_type

Should write to results/tables/:
    <project>_ranked_genes.csv   gene_id, gene_name, log2FC, p_value, padj, status

TODO (Daniel): log2 fold change, per-gene t-test, multiple-testing
correction, ranked table, volcano plots.
"""

from pathlib import Path

PROCESSED_DIR = Path("data") / "processed"
TABLES_DIR = Path("results") / "tables"
FIGURES_DIR = Path("results") / "figures"
