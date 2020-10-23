
# Description -------------------------------------------------------------

# Run GSVA using tidy-signatures.Rds and IMvigor210


# Prepare Workspace -------------------------------------------------------

library(readr)
library(DESeq2)
library(GSVA)


# Load Dataframes ---------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/norm_dds.Rds")
signatures <- read_rds("./data/imvigor210/sig-with-null.Rds")


# Run GSVA ----------------------------------------------------------------

per_tumor_signature <- gsva(assay(imvigor, 2), signatures, mx.diff = T, kcdf = "Poisson")


# Export GSVA -------------------------------------------------------------

write_rds(per_tumor_signature, "./data/imvigor210/gsva-scores.Rds")
