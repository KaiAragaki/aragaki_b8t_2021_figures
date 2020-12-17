# Tidy signatures for GSVA:

# Pan B-Cell Signature
# Source: https://jitc.bmj.com/content/5/1/18.long

# CD8+ T Effector Signature (cd8_rosen): 
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5480242/#SD1
# Supplementary figure 6

# CD8+ T (cd8_prat):
# Source: https://cancerres.aacrjournals.org/content/canres/77/13/3540/F1.large.jpg?width=800&height=600&carousel=1

# CD8+ T (cd8_sade):
# Source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6641984/bin/NIHMS1510803-supplement-11.xlsx
# CD8_G

# CD8+ T (cd8_fehr):
# Source: https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(16)00587-0/fulltext

# IFNg Signature
# Source: https://www.jci.org/articles/view/91190/table/2


# Prepare Workspace -------------------------------------------------------

library(dplyr)
library(readxl)


# Load Signatures ---------------------------------------------------------

b_cell <- c("BLK", "CD19", "FCRL2", "MS4A1", "KIAA0125", "TNFRSF17", 
            "TCL1A", "SPIB", "PNOC")

cd8_rosen <- c("CD8A", "GZMA", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1", 
               "TBX21")

cd8_prat <- c("PRF1", "CD8A", "CD8B", "GZMM", "FLT3LG")

cd8_sade <- c("IL7R", "GPR183", "LMNA", "NR4A3", "TCF7", "MGAT4A", "CD55", 
              "AIM1", "PER1", "FOSL2", "EGFR1", "TSPYL2", "YPEL5", 
              "CSRNP1", "REL", "SKIL", "PIK3R1", "FOXP1", "RGCC", "PFKFB3",
              "MYADM", "ZFP36L2", "USP36", "TC2N", "FAM177A1", "BTG2", 
              "TSC22D2", "FAM65B", "STAT4", "RGP5", "NEU1", "IRFD1", 
              "PDE4B", "NR4A1")

cd8_fehr <- c("CD8A", "EOMES", "PRF1", "IFNG", "CD274")

ifng <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")


# Put into List -----------------------------------------------------------

signatures <- list(b_cell = b_cell, cd8_rose = cd8_rosen, 
                   cd8_prat = cd8_prat, cd8_fehr = cd8_fehr, ifng = ifng)

# Write -------------------------------------------------------------------

write_rds(signatures, "./data/signatures/signatures.Rds")
