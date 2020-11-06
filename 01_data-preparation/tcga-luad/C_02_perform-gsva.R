
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA SKCM


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

skcm <- read_rds("./data/tcga-luad/A_02_normalized-counts.Rds")
signatures <- read_rds("./data/tcga-luad/C_01_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

per_tumor_sig <-  gsva(assay(skcm, 2), signatures, mx.diff = T, kcdf = "Gaussian")


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_sig, "./data/tcga-luad/C_02_gsva-scores.Rds")
