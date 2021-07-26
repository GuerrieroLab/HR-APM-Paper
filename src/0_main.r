setwd("~/Dropbox (HMS)/_BTIL_/progress/hr_apm_paper/src/")

## 1. create gene signatures from the excel files
source("gene_list_prep.r")

## 2. Prepare TCGA data
source("tcga_prep.r")

## 3. Prepare METABRIC data
source("metabric_prep.r")

## 4. Correlations in TCGA and METABRIC
source("tcga_metabric_correlation.r")

## 5. Looking into the overlap between the clusters
source("overlap_tcga_metabric.r")

## 6. Looking into the correlation in different cancer types
source("bc_subtype_panel_genes.r")



