
# Description -------------------------------------------------------------

# Run GSVA on Riaz


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

riaz <- read_rds("./data/riaz/dds.Rds")
sigs <- read_rds("./data/riaz/02_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

pre <- riaz[, riaz$pre_or_on_treat == "Pre"]

on  <- riaz[, riaz$pre_or_on_treat == "On"]

riaz_sig <- gsva(assay(riaz, 2), sigs, mx.diff = T, kcdf = "Gaussian")
pre_sig <-  gsva(assay(pre, 2), sigs, mx.diff = T, kcdf = "Gaussian")
on_sig <-   gsva(assay(on, 2), sigs, mx.diff = T, kcdf = "Gaussian")

spec_sigs <- cbind(pre_sig, on_sig)
rownames(spec_sigs) <- paste0(rownames(spec_sigs), "_spec")
riaz_sig <- riaz_sig[,match(colnames(spec_sigs), colnames(riaz_sig))]
per_tumor_sig <- rbind(spec_sigs, riaz_sig)

# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_sig, "./data/riaz/04_gsva-scores.Rds")
