
# Description -------------------------------------------------------------

# Matches GSVA scores to their respective tumors in the colData


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(DESeq2)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/norm_dds.Rds")

imvigor_colData <- imvigor %>%
        colData() %>%
        as_tibble(rownames = "sample")

gsva <- read_rds("./data/imvigor210/gsva-scores.Rds") %>%
        t() %>%
        as_tibble(rownames = "sample")

joined <- inner_join(imvigor_colData, gsva, by = "sample")

colData(imvigor) <- DataFrame(joined)

# Write -------------------------------------------------------------------

write_rds(imvigor, "./data/imvigor210/dds-gsva.Rds")
