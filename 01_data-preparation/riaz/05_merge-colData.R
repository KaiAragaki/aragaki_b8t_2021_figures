
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)


# Read in Data ------------------------------------------------------------

riaz <- read_rds("./data/riaz/dds.Rds")

riaz_col_data <- riaz %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/riaz/04_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(riaz_col_data, gsva, by = "sample")

joined <- as.data.frame(joined, row.names = riaz_col_data$sample)


# Write -------------------------------------------------------------------

write_rds(joined, "./data/riaz/05_merge-colData.rds")
