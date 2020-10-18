
# Description -------------------------------------------------------------

# Making GSVA-compatible signatures from MSigDB Hallmark gene sets
# gene sets can be found at https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=H


# Prepare Workspace -------------------------------------------------------

library(tidyverse)


# Define Functions --------------------------------------------------------

# Determine if header of file has '>', denoting MSigDB file
check_msigdb <- function(path) {
        read_tsv(file = path, n_max = 1, col_types = 'c') %>%
                unlist() %>%
                str_detect("^>")
}

read_no_description <- function(path) {
        read_tsv(path, col_types = "c") %>%
                dplyr::slice(-1) %>%
                pull()
}


# Scan through folder -----------------------------------------------------

# Makes checks to ensure msigdb file
sets <- list.files("./data/signatures") %>%
        as_tibble() %>%
        filter(str_detect(value, ".txt$")) %>%
        mutate(path = paste0("./data/signatures/", value)) %>%
        mutate(msig = map_lgl(path, check_msigdb),
               name = str_remove(value, ".txt$")) %>%
        filter(msig)


# Read in Signatures ------------------------------------------------------

signatures <- map(sets$path, read_no_description)
names(signatures) <- sets$name


# Write Signature File ----------------------------------------------------

write_rds(signatures, "./data/signatures/msig-signatures.Rds")
