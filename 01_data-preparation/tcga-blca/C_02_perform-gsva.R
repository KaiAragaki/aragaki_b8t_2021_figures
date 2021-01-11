
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA BLCA


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

blca <- read_rds("./data/tcga-blca/A_02_normalized-counts.Rds")
signatures <- read_rds("./data/tcga-blca/C_02_signatures-with-null.Rds")


# Ensure Genes are in Dataset ----------------------------------------

signatures$b_cell[5] <- "FAM30A"

genes <- unlist(signatures) %>% 
        as_tibble()



if(nrow(genes[which(!(genes$value %in% rownames(blca))),]) == 0) {
        print("All signature genes have match in dataset.")
}


# Run GSVA ----------------------------------------------------------------

per_tumor_signature <- gsva(assay(blca,2), signatures, mx.diff = T, 
                            kcdf = "Gaussian")


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_signature, "./data/tcga-blca/C_03_gsva-scores.Rds")
