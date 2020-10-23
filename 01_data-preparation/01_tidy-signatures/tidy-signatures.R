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

sig_b_cell <- read_xlsx("./data/signatures/DanaherSig.xlsx") %>%
        dplyr::pull(1)

sig_cd8_t_eff <- c("CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", "TBX21")


# Put into List -----------------------------------------------------------

signatures <- list(b_cell = sig_b_cell, cd8_t_eff = sig_cd8_t_eff)

# Write -------------------------------------------------------------------

write_rds(signatures, "./data/signatures/signatures.Rds")
