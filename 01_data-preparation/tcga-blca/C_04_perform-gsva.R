
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA BLCA


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

blca <- read_rds("./data/tcga-blca/A_02_normalized-counts.Rds")
signatures <- read_rds("./data/tcga-blca/C_03_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

per_tumor_signature <- gsva(assay(blca,2), signatures, mx.diff = T, kcdf = "Poisson")


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_signature, "./data/tcga-blca/C_04_gsva-scores.Rds")
