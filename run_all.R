# Run from the project root: Rscript run_all.R
# Downloads and prepares the TCGA data, then exports it for python/dgea.py.

if (!file.exists(file.path("R", "load_data.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "load_data.R"))

PROJECTS <- c("TCGA-BRCA", "TCGA-LUAD")

for (project in PROJECTS) {
  export_for_python(load_cancer(project))
}

message("Data ready in ", PROCESSED_DIR, ". Next: python python/dgea.py")
