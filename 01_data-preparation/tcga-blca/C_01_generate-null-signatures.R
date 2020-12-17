
# Description -------------------------------------------------------------

# Generate Null Signatures for GSVA


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)


# Define Functions --------------------------------------------------------

generate_null_sig <- function(names, signature) {
        sig_name <- enexpr(names)
        # Allows signature null distributions to be different between
        # signatures, but the same between runs
        sig_seed <- substitute(sig_name) %>%
                deparse() %>%
                charToRaw() %>%
                as.numeric() %>%
                sum()
        sample(rownames(blca), size = length(signature))
}


# Load Data ---------------------------------------------------------------

blca <- read_rds("./data/tcga-blca/A_02_normalized-counts.Rds")


# Remove Blank HGNC IDs ---------------------------------------------------

blca <- blca[(rownames(blca) != "" & !is.na(rownames(blca))), ]


# Load Signatures ---------------------------------------------------------

signatures <- read_rds("./data/signatures/signatures.Rds")


# Generate Null Signatures ------------------------------------------------

signatures_null <- map2(names(signatures), signatures, generate_null_sig)
names(signatures_null) <- paste0(names(signatures), "_null")
signatures <- c(signatures, signatures_null)


# Export ------------------------------------------------------------------

write_rds(signatures, "./data/tcga-blca/C_02_signatures-with-null.Rds")
