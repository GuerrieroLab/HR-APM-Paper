## Relationship between ER expression and APM signature

This repository contains code to explore anticorrelations between the ESR and APM/TC modules (Cui, Shimada, Goldbertg *et al.*)

1. create gene signatures from the excel files - [1_gene_list_prep.Rmd](src/1_gene_list_prep.Rmd)
2. Prepare TCGA data - [2_tcga_prep.Rmd](src/2_tcga_prep.Rmd)
3. Prepare METABRIC data - [3_metabric_prep.Rmd](src/3_metabric_prep.Rmd)
4. Correlations in TCGA and METABRIC - [tcga_metabric_correlation.Rmd](src/tcga_metabric_correlation.Rmd)
5. Looking into the overlap between the clusters - [overlap_tcga_metabric.Rmd](src/overlap_tcga_metabric.Rmd)
