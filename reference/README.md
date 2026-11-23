# Reference data

`known_cancer_genes.csv`: 638 genes with clinical evidence in cancer, taken from
the gene summaries of CIViC (Clinical Interpretation of Variants in Cancer,
https://civicdb.org), downloaded 2026-11-23 from
https://civicdb.org/downloads/nightly/nightly-GeneSummaries.tsv, keeping rows of
type `Gene` and the columns `name` (as `gene_name`) and `entrez_id`.

CIViC data is released under CC0 1.0 (public domain), so it can be redistributed
here. Griffith et al., Nature Genetics 2017, doi:10.1038/ng.3774.

Used by `R/validate.R`.
