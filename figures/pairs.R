
# Description -------------------------------------------------------------

# Compare Signature Scores


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(DESeq2)
library(viridisLite)
library(GGally)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")

gsva_scores <- imvigor %>% 
        colData(imvigor) %>% 
        as_tibble() %>% 
        select(b_cell:ifng)

ggpairs(gsva_scores)

