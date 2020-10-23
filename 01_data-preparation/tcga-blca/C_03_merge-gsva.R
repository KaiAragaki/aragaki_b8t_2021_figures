
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)


# Read in Data ------------------------------------------------------------

blca <- read_rds("./data/tcga-blca/A_02_normalized-counts.Rds")

blca_coldata <- blca %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/tcga-blca/C_03_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(blca_coldata, gsva, by = "sample")

colData(blca) <- DataFrame(joined)

# Write -------------------------------------------------------------------

write_rds(blca, "./data/tcga-blca/C_04_merge-gsva.rds")
