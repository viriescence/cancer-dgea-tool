# Run from the project root.
#   Rscript run_all.R         uses the committed sample in data/sample/
#   Rscript run_all.R --full  downloads full TCGA data (~10 GB), rebuilds the sample
# Either way, writes the CSV files that python/dgea.py reads to data/processed/.

if (!file.exists(file.path("R", "load_data.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "load_data.R"))

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
}

message("Data ready in ", PROCESSED_DIR, ". Next: python python/dgea.py")
