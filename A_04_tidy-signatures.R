# Tidy signatures for GSVA:

# Pan B-Cell Signature
# Danaher et al.

# CD8+ T_eff
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480242/#SD1
# Supplementary figure 6


# Prepare Workspace -------------------------------------------------------

library(dplyr)
library(readxl)


# Load Signatures ---------------------------------------------------------

#sig_b_cell <- read_xlsx("./data/signatures/DanaherSig.xlsx") %>%
 # dplyr::pull(1)

sig_b_cell <- c("BLK","CD19","FCRL2","MS4A1","KIAA0125","TNFRSF17","TCL1A","SPIB","PNOC")

sig_cd8_t_eff <- c("CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21")


# Put into List -----------------------------------------------------------

signatures <- list(b_cell = sig_b_cell, cd8_t_eff = sig_cd8_t_eff)
View(signatures)
# Write -------------------------------------------------------------------

write_rds(signatures, "H:/Rusty_B8T/data/TCGA_BLCA/signatures/signatures.Rds")