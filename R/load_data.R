# Download, label, filter and normalize TCGA expression data.

suppressPackageStartupMessages(library(edgeR))

RAW_DIR       <- file.path("data", "raw")
PROCESSED_DIR <- file.path("data", "processed")
SAMPLE_DIR    <- file.path("data", "sample")

# Only the full TCGA download needs these packages; the committed sample does not.
require_bioc <- function() {
  for (pkg in c("TCGAbiolinks", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required for the full TCGA download: ",
           "BiocManager::install(\"", pkg, "\")")
    }
  }
}

download_tcga <- function(project, raw_dir = RAW_DIR) {
  require_bioc()
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

  query <- TCGAbiolinks::GDCquery(
    project       = project,
    data.category = "Transcriptome Profiling",
    data.type     = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  TCGAbiolinks::GDCdownload(query, directory = raw_dir)
  TCGAbiolinks::GDCprepare(query, directory = raw_dir)
}

label_samples <- function(se) {
  meta <- as.data.frame(SummarizedExperiment::colData(se))

  group <- rep(NA_character_, nrow(meta))
  group[meta$sample_type == "Primary Tumor"]       <- "tumor"
  group[meta$sample_type == "Solid Tissue Normal"] <- "normal"

  # One sample per patient per group, so samples are independent.
  patient <- substr(colnames(se), 1, 12)
  keep <- !is.na(group) & !duplicated(paste(patient, group))

  se <- se[, keep]
  se$group <- factor(group[keep], levels = c("normal", "tumor"))
  se
}

extract_counts <- function(se, protein_coding_only = TRUE) {
  counts <- SummarizedExperiment::assay(se, "unstranded")
  info   <- as.data.frame(SummarizedExperiment::rowData(se))[, c("gene_id", "gene_name", "gene_type")]
  rownames(info) <- info$gene_id

  if (protein_coding_only) {
    pc <- info$gene_type == "protein_coding"
    counts <- counts[pc, , drop = FALSE]
    info   <- info[pc, , drop = FALSE]
  }
  list(counts = counts, genes = info)
}

# Keep genes with CPM > min_cpm in at least as many samples as the smallest group.
filter_low_expression <- function(counts, group, min_cpm = 1) {
  min_samples <- min(table(group))
  keep <- rowSums(cpm(counts) > min_cpm) >= min_samples
  counts[keep, , drop = FALSE]
}

# TMM-normalized log2 CPM.
normalize_log_cpm <- function(counts, group) {
  dge <- DGEList(counts = counts, group = group)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 1)
}

# Full TCGA download. Returns list(project, counts, logcpm, group, genes);
# cached in data/processed/.
load_cancer <- function(project, force = FALSE) {
  cache <- file.path(PROCESSED_DIR, paste0(project, ".rds"))
  if (file.exists(cache) && !force) {
    message("Using cached ", cache)
    return(readRDS(cache))
  }

  message("[", project, "] downloading / loading raw counts ...")
  se <- label_samples(download_tcga(project))

  parts  <- extract_counts(se)
  group  <- se$group
  counts <- filter_low_expression(parts$counts, group)
  genes  <- parts$genes[rownames(counts), ]

  stopifnot(identical(colnames(counts), colnames(se)))
  if (min(table(group)) < 3) {
    stop(project, " has fewer than 3 samples in one group (",
         paste(names(table(group)), table(group), collapse = ", "),
         "); a tumor vs. normal comparison is not possible.")
  }

  message("[", project, "] ", table(group)[["tumor"]], " tumor, ",
          table(group)[["normal"]], " normal samples; ",
          nrow(counts), " genes after filtering")

  out <- list(
    project = project,
    counts  = counts,
    logcpm  = normalize_log_cpm(counts, group),
    group   = group,
    genes   = genes
  )

  dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out, cache)
  out
}

# All normals plus a seeded random subset of tumors, small enough to commit.
make_sample <- function(data, n_tumor = 150, seed = 123) {
  set.seed(seed)
  tumors <- which(data$group == "tumor")
  keep <- sort(c(which(data$group == "normal"),
                 sample(tumors, min(n_tumor, length(tumors)))))

  list(project = data$project,
       counts  = data$counts[, keep],
       group   = data$group[keep],
       genes   = data$genes)
}

# Writes <project>_counts.csv.gz, <project>_groups.csv, <project>_genes.csv.
write_sample <- function(sample, dir = SAMPLE_DIR) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  p <- sample$project

  con <- gzfile(file.path(dir, paste0(p, "_counts.csv.gz")), "w")
  write.csv(sample$counts, con, quote = FALSE)
  close(con)

  write.csv(data.frame(sample_id = colnames(sample$counts),
                       group     = as.character(sample$group)),
            file.path(dir, paste0(p, "_groups.csv")),
            row.names = FALSE, quote = FALSE)

  write.csv(sample$genes[, c("gene_id", "gene_name", "gene_type")],
            file.path(dir, paste0(p, "_genes.csv")),
            row.names = FALSE)
}

# Same return shape as load_cancer(), built from the committed sample.
load_sample <- function(project, dir = SAMPLE_DIR) {
  counts <- as.matrix(read.csv(file.path(dir, paste0(project, "_counts.csv.gz")),
                               row.names = 1, check.names = FALSE))
  groups <- read.csv(file.path(dir, paste0(project, "_groups.csv")))
  genes  <- read.csv(file.path(dir, paste0(project, "_genes.csv")),
                     row.names = NULL)
  if (!"gene_name" %in% names(genes)) genes$gene_name <- genes$gene_id
  if (!"gene_type" %in% names(genes)) genes$gene_type <- "unknown"
  rownames(genes) <- genes$gene_id

  group <- factor(groups$group[match(colnames(counts), groups$sample_id)],
                  levels = c("normal", "tumor"))
  if (anyNA(group)) {
    stop("Every sample in ", project, "_counts.csv.gz needs a 'tumor' or ",
         "'normal' label in ", project, "_groups.csv.")
  }
  stopifnot(setequal(rownames(counts), genes$gene_id))

  list(project = project,
       counts  = counts,
       logcpm  = normalize_log_cpm(counts, group),
       group   = group,
       genes   = genes[rownames(counts), ])
}

# Plain-text copy of a loaded dataset for python/dgea.py:
# <project>_logcpm.csv.gz, <project>_groups.csv, <project>_genes.csv
export_for_python <- function(data, out_dir = PROCESSED_DIR) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p <- data$project

  con <- gzfile(file.path(out_dir, paste0(p, "_logcpm.csv.gz")), "w")
  write.csv(round(data$logcpm, 4), con, quote = FALSE)
  close(con)

  write.csv(data.frame(sample_id = colnames(data$logcpm),
                       group     = as.character(data$group)),
            file.path(out_dir, paste0(p, "_groups.csv")),
            row.names = FALSE, quote = FALSE)

  write.csv(data$genes[, c("gene_id", "gene_name", "gene_type")],
            file.path(out_dir, paste0(p, "_genes.csv")),
            row.names = FALSE)

  invisible(file.path(out_dir, paste0(p, "_*")))
}
