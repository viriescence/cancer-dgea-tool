# DESeq2 differential expression (tumor vs. normal) on raw counts.

run_deseq2 <- function(data, lfc_cutoff = 1, fdr_cutoff = 0.05) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required: BiocManager::install(\"DESeq2\")")
  }

  counts <- data$counts
  storage.mode(counts) <- "integer"
  coldata <- data.frame(group = data$group, row.names = colnames(counts))

  dds <- DESeq2::DESeqDataSetFromMatrix(counts, coldata, design = ~ group)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", "tumor", "normal"))

  out <- data.frame(
    gene_id   = rownames(res),
    gene_name = data$genes[rownames(res), "gene_name"],
    log2FC    = res$log2FoldChange,
    p_value   = res$pvalue,
    padj      = res$padj,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  significant <- !is.na(out$padj) & out$padj < fdr_cutoff
  out$status <- ifelse(significant & out$log2FC >=  lfc_cutoff, "up",
                ifelse(significant & out$log2FC <= -lfc_cutoff, "down", "ns"))

  out[order(out$padj, -abs(out$log2FC)), ]
}
