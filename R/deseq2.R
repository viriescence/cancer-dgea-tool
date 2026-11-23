# DESeq2 differential expression (tumor vs. normal) on raw counts.

run_deseq2 <- function(data, lfc_cutoff = 1, fdr_cutoff = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required: BiocManager::install(\"DESeq2\")")
  }

  # DESeq2 wants raw integer counts and a table describing each sample.
  counts <- data$counts
  storage.mode(counts) <- "integer"
  coldata <- data.frame(group = data$group, row.names = colnames(counts))

  # Fit the model (counts ~ group) and extract tumor vs. normal.
  # Normal is the first factor level, so a positive log2FC means higher in tumor.
  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", "tumor", "normal"))

  # Same columns as the Python table, so the two can be compared directly.
  out <- data.frame(
    gene_id   = rownames(res),
    gene_name = data$genes[rownames(res), "gene_name"],
    log2FC    = res$log2FoldChange,
    p_value   = res$pvalue,
    padj      = res$padj,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # up/down need both significance (padj) and a big enough change (log2FC).
  significant <- !is.na(out$padj) & out$padj < fdr_cutoff
  out$status <- ifelse(significant & out$log2FC >=  lfc_cutoff, "up",
                ifelse(significant & out$log2FC <= -lfc_cutoff, "down", "ns"))

  # Most significant genes first; ties broken by size of the change.
  out[order(out$padj, -abs(out$log2FC)), ]
}
