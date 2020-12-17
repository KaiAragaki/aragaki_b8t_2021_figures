
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
  sample(rownames(blca), size = length(signature))
}


# Load Data ---------------------------------------------------------------

blca <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/A_02_normalized-counts.Rds")



# Remove Blank HGNC IDs ---------------------------------------------------

blca <- blca[(rownames(blca) != "" & !is.na(rownames(blca))), ]


# Load Signatures ---------------------------------------------------------

signatures <- read_rds("H:/Rusty_B8T/data/TCGA_BLCA/signatures/signatures.Rds")


# Generate Null Signatures ------------------------------------------------
names(signatures)
signatures_null <- map2(names(signatures), signatures, generate_null_sig)
?map2
names(signatures_null) <- paste0(names(signatures), "_null")
signatures <- c(signatures, signatures_null)

View(signatures_null)
View(signatures)

# Export ------------------------------------------------------------------

write_rds(signatures, "H:/Rusty_B8T/data/TCGA_BLCA/C_02_signatures-with-null.Rds")