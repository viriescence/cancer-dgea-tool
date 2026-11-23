# Run from the project root.
#   Rscript run_all.R                  BRCA and LUAD from the committed sample
#   Rscript run_all.R TCGA-COAD ...    any TCGA project (or your own dataset in
#                                      data/sample/, see README); downloads and
#                                      samples it first if no sample exists
#   Rscript run_all.R --full           rebuild the sample(s) from TCGA (~10 GB each)
# Writes the CSV files that python/dgea.py reads to data/processed/, plus the
# DESeq2 table and volcano plot for each project to results/.

if (!file.exists(file.path("R", "load_data.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "load_data.R"))
source(file.path("R", "deseq2.R"))
source(file.path("R", "plots.R"))

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args
projects <- setdiff(args, "--full")
if (length(projects) == 0) projects <- c("TCGA-BRCA", "TCGA-LUAD")

for (project in projects) {
  has_sample <- file.exists(file.path(SAMPLE_DIR, paste0(project, "_counts.csv.gz")))
  if (full || !has_sample) {
    write_sample(make_sample(load_cancer(project)))
  }
  data <- load_sample(project)
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
