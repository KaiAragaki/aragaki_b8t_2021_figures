
# Description -------------------------------------------------------------

# Run GSVA on Liu dataset


# Prepare Workspace -------------------------------------------------------

library(readr)
library(SummarizedExperiment)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

liu <- read_rds("./data/liu/se.Rds")
sigs <- read_rds("./data/liu/02_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

ctla4_prior <- liu[, liu$priorCTLA4 == 1]
ctla4_naive <- liu[, liu$priorCTLA4 == 0]

liu_sigs        <-  gsva(assay(liu, 2), sigs, mx.diff = T, kcdf = "Gaussian")
ctla4_prior_sig <-  gsva(assay(ctla4_prior, 2), sigs, mx.diff = T, kcdf = "Gaussian")
ctla4_naive_sig <-  gsva(assay(ctla4_naive, 2), sigs, mx.diff = T, kcdf = "Gaussian")

spec_sigs <- cbind(ctla4_prior_sig, ctla4_naive_sig)
rownames(spec_sigs) <- paste0(rownames(spec_sigs), "_spec")
liu_sigs <- liu_sigs[,match(colnames(spec_sigs), colnames(liu_sigs))]
per_tumor_sig <- rbind(spec_sigs, liu_sigs)


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_sig, "./data/liu/03_gsva-scores.Rds")
