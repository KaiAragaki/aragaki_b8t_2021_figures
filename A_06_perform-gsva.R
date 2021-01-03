
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA BLCA


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

blca <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/A_02_normalized-counts.Rds")
signatures <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/C_02_signatures-with-null.Rds")

# Ensure Genes are in Dataset ----------------------------------------

genes <- unlist(signatures) %>% 
  as_tibble()


if(nrow(genes[which(!(genes$value %in% rownames(blca))),]) == 0) {
  print("All signature genes have match in dataset.")
}


# Run GSVA changed possion to Gaussian----------------------------------------------------------------

per_tumor_signature <- gsva(assay(blca,2), signatures, mx.diff = T, kcdf = "Gaussian")

View(per_tumor_signature)
# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_signature, "H:/Rusty_B8T/data/TCGA_BLCA/C_03_gsva-scores.Rds")