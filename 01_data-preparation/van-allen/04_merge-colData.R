
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)


# Read in Data ------------------------------------------------------------

va <- read_rds("./data/van-allen/se.rds")

va_col_data <- va %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/van-allen/03_gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(va_col_data, gsva, by = "sample")

joined <- as.data.frame(joined, row.names = va_col_data$sample)


# Write -------------------------------------------------------------------

write_rds(joined, "./data/van-allen/04_merge-colData.rds")
