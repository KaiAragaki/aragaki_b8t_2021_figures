
# Description -------------------------------------------------------------

# Run GSVA on van Allen dataset


# Prepare Workspace -------------------------------------------------------

library(readr)
library(SummarizedExperiment)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

va <- read_rds("./data/van-allen/se.Rds")
sigs <- read_rds("./data/van-allen/02_signatures-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

scores <-  gsva(assay(va, 2), sigs, mx.diff = T, kcdf = "Gaussian")


# Export GSVA -------------------------------------------------------------

write_rds(scores, "./data/van-allen/03_gsva-scores.Rds")
