
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA BLCA


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

blca <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/A_02_normalized-counts.Rds")
signatures <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/C_02_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

per_tumor_signature <- gsva(assay(blca,2), signatures, mx.diff = T, kcdf = "Poisson")

View(per_tumor_signature)
# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_signature, "H:/Rusty_B8T/data/TCGA_BLCA/C_03_gsva-scores.Rds")