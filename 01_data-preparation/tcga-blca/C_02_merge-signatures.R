
# Description -------------------------------------------------------------

# Merges Standard and Nonstandard Signatures


# Prepare Workspace -------------------------------------------------------

library(readr)


# Read in Datasets --------------------------------------------------------

msig <- read_rds("./data/signatures/msig-signatures.Rds")
nonstandard <- read_rds("./data/signatures/nonstandard-signatures.Rds")


# Join --------------------------------------------------------------------

signatures <- c(msig, nonstandard)


# Write -------------------------------------------------------------------

write_rds(signatures, "./data/signatures/merged-signatures.Rds")
