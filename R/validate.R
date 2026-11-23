# Validation: are known cancer genes (reference/known_cancer_genes.csv) enriched
# among the genes each method calls significant?

KNOWN_GENES_FILE <- file.path("reference", "known_cancer_genes.csv")

load_known_genes <- function(path = KNOWN_GENES_FILE) {
  unique(read.csv(path, stringsAsFactors = FALSE)$gene_name)
}

# res: ranked table with gene_name, log2FC, padj, status (DESeq2 or Python).
# Two questions, each tested with Fisher's exact test:
#   1. Are known genes more often significant than other genes?
#   2. Are known genes over-represented among the top_n most significant genes?
validate_known_genes <- function(res, known, top_n = 500) {
  res <- res[!is.na(res$padj), ]
  res <- res[order(res$padj, -abs(res$log2FC)), ]

  is_known <- res$gene_name %in% known
  significant <- res$status != "ns"
  in_top <- seq_len(nrow(res)) <= top_n

  data.frame(
    known_tested     = sum(is_known),
    known_significant = sum(is_known & significant),
    pct_known_significant = round(100 * mean(significant[is_known]), 1),
    pct_other_significant = round(100 * mean(significant[!is_known]), 1),
    p_significant    = fisher.test(table(is_known, significant))$p.value,
    known_in_top     = sum(is_known & in_top),
    # How many known genes we would expect in the top list by chance alone.
    expected_in_top  = round(top_n * mean(is_known), 1),
    p_top            = fisher.test(table(is_known, in_top))$p.value
  )
}

# Runs on every ranked table that exists in results/tables/.
# Missing tables (for example before the Python step has run) are skipped.
validate_all <- function(projects = c("TCGA-BRCA", "TCGA-LUAD"), top_n = 500) {
  known <- load_known_genes()
  methods <- c(DESeq2 = "_deseq2.csv", Python = "_ranked_genes.csv")

  rows <- list()
  for (project in projects) {
    for (method in names(methods)) {
      path <- file.path("results", "tables", paste0(project, methods[[method]]))
      if (!file.exists(path)) next
      out <- validate_known_genes(read.csv(path, stringsAsFactors = FALSE), known, top_n)
      rows[[length(rows) + 1]] <- cbind(project = project, method = method, out)
    }
  }

  # One row per project and method; rounded so the CSV stays readable.
  summary <- do.call(rbind, rows)
  summary$p_significant <- signif(summary$p_significant, 3)
  summary$p_top <- signif(summary$p_top, 3)

  dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, file.path("results", "tables", "known_gene_validation.csv"),
            row.names = FALSE)
  summary
}
