
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and TCGA SKCM


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

skcm <- read_rds("./data/tcga-skcm/A_02_normalized-counts.Rds")
signatures <- read_rds("./data/tcga-skcm/C_01_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

na         <- skcm[,  is.na(skcm$patient.tumor_tissue_site)]
met        <- skcm[, !is.na(skcm$patient.tumor_tissue_site) & skcm$patient.tumor_tissue_site == "distant metastasis"]
prime      <- skcm[, !is.na(skcm$patient.tumor_tissue_site) & skcm$patient.tumor_tissue_site == "primary tumor"]
reg_tissue <- skcm[, !is.na(skcm$patient.tumor_tissue_site) & skcm$patient.tumor_tissue_site == "regional cutaneous or subcutaneous tissue (includes satellite and in-transit metastasis)"]
lymph      <- skcm[, !is.na(skcm$patient.tumor_tissue_site) & skcm$patient.tumor_tissue_site == "regional lymph node"]

per_tumor_sig  <- gsva(assay(skcm, 2), signatures, mx.diff = T, kcdf = "Gaussian")
na_sig         <- gsva(assay(na, 2), signatures, mx.diff = T, kcdf = "Gaussian")
met_sig        <- gsva(assay(met, 2), signatures, mx.diff = T, kcdf = "Gaussian")
prime_sig      <- gsva(assay(prime, 2), signatures, mx.diff = T, kcdf = "Gaussian")
reg_tissue_sig <- gsva(assay(reg_tissue, 2), signatures, mx.diff = T, kcdf = "Gaussian")
lymph_sig      <- gsva(assay(lymph, 2), signatures, mx.diff = T, kcdf = "Gaussian")

spec_sigs <- cbind(na_sig, met_sig, prime_sig, reg_tissue_sig, lymph_sig)
rownames(spec_sigs) <- paste0(rownames(spec_sigs), "_spec")
per_tumor_sig <- per_tumor_sig[,match(colnames(spec_sigs), colnames(per_tumor_sig))]
per_tumor_sig <- rbind(spec_sigs, per_tumor_sig)


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_sig, "./data/tcga-skcm/C_02_gsva-scores.Rds")
