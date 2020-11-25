
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(DESeq2)


# Read in Data ------------------------------------------------------------


luad <- read_rds("./data/tcga-luad/A_02_normalized-counts.Rds")

luad_coldata <- luad %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/tcga-luad/C_02_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(skcm_coldata, gsva, by = "sample")

joined <- as.data.frame(joined, row.names = skcm_coldata$sample)


# Write -------------------------------------------------------------------

write_rds(joined, "./data/tcga-luad/C_03_merge-gsva.rds")
