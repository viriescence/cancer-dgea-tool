# Run from the project root, after run_all.R and python/dgea.py.
# Checks both methods' results against reference/known_cancer_genes.csv and
# writes results/tables/known_gene_validation.csv.

if (!file.exists(file.path("R", "validate.R"))) {
  stop("Run this script from the project root.")
}

source(file.path("R", "validate.R"))
print(validate_all(), row.names = FALSE)
