# Run from the project root.
#   Rscript run_all.R         uses the committed sample in data/sample/
#   Rscript run_all.R --full  downloads full TCGA data (~10 GB), rebuilds the sample
# Writes the CSV files that python/dgea.py reads to data/processed/, plus the
# DESeq2 table and volcano plot for each cancer to results/.

if (!file.exists(file.path("R", "load_data.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "load_data.R"))
source(file.path("R", "deseq2.R"))
source(file.path("R", "plots.R"))

PROJECTS <- c("TCGA-BRCA", "TCGA-LUAD")
full <- "--full" %in% commandArgs(trailingOnly = TRUE)

for (project in PROJECTS) {
  if (full) {
    data <- load_cancer(project)
    write_sample(make_sample(data))
  } else {
    data <- load_sample(project)
  }
  export_for_python(data)

  res <- run_deseq2(data)
  dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
  write.csv(res, file.path("results", "tables", paste0(project, "_deseq2.csv")),
            row.names = FALSE)
  save_volcano(res,
               file.path("results", "figures", paste0(project, "_deseq2_volcano.png")),
               paste(project, "- DESeq2"))
}

message("Done. Python input in ", PROCESSED_DIR, "; DESeq2 results in results/.")
