# Volcano plot of a ranked gene table (columns gene_name, log2FC, padj, status).

save_volcano <- function(res, path, title, top_n = 10,
                         lfc_cutoff = 1, fdr_cutoff = 0.05) {
  res <- res[!is.na(res$padj), ]
  res$neg_log10_padj <- -log10(pmax(res$padj, 1e-300))

  labelled <- head(res[res$status != "ns", ], top_n * 2)

  p <- ggplot2::ggplot(res, ggplot2::aes(log2FC, neg_log10_padj, colour = status)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    ggplot2::geom_text(data = labelled, ggplot2::aes(label = gene_name),
                       size = 2.8, vjust = -0.6, check_overlap = TRUE,
                       show.legend = FALSE) +
    ggplot2::scale_colour_manual(
      values = c(up = "#c0392b", down = "#2471a3", ns = "grey70")) +
    ggplot2::labs(title = title, x = "log2 fold change (tumor vs. normal)",
                  y = "-log10 adjusted p-value", colour = NULL) +
    ggplot2::theme_bw()

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, p, width = 7, height = 6, dpi = 150)
}
