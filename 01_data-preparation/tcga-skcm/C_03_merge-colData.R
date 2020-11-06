
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)


# Read in Data ------------------------------------------------------------

skcm <- read_rds("./data/tcga-skcm/A_02_normalized-counts.Rds")

skcm_coldata <- skcm %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/tcga-skcm/C_02_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(skcm_coldata, gsva, by = "sample")

joined <- as.data.frame(joined, row.names = skcm_coldata$sample)


# Write -------------------------------------------------------------------

write_rds(joined, "./data/tcga-skcm/C_03_merge-colData.rds")
