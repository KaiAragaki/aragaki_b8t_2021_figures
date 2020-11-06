
# Description -------------------------------------------------------------

# Generate Null Signatures for GSVA


# Prepare Workspace -------------------------------------------------------

library(tidyverse)


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
        sample(rownames(liu), size = length(signature))
}


# Load Data ---------------------------------------------------------------

liu <- read_rds("./data/liu/se.Rds")


# Remove Blank HGNC IDs ---------------------------------------------------

liu <- liu[(rownames(liu) != "" & !is.na(rownames(liu))), ]


# Load Signatures ---------------------------------------------------------

signatures <- read_rds("./data/signatures/signatures.Rds")


# Generate Null Signatures ------------------------------------------------

signatures_null <- map2(names(signatures), signatures, generate_null_sig)
names(signatures_null) <- paste0(names(signatures), "_null")
signatures <- c(signatures, signatures_null)


# Export ------------------------------------------------------------------

write_rds(signatures, "./data/liu/02_signatures-with-null.Rds")
