
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)


# Read in Data ------------------------------------------------------------

liu <- read_rds("./data/liu/se.Rds")

col_data <- liu %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/liu/03_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(col_data, gsva, by = "sample")

joined <- as.data.frame(joined, row.names = col_data$sample)


# Write -------------------------------------------------------------------

write_rds(joined, "./data/liu/04_merge-colData.Rds")
