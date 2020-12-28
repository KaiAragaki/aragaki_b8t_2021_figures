

# Description -------------------------------------------------------------

# Scrape FFPE cases from Broad Firehose


# Prepare Workspace -------------------------------------------------------

library(rvest)
library(readr)
library(dplyr)


# Read in Data ------------------------------------------------------------

ffpe <- read_html("http://gdac.broadinstitute.org/runs/stddata__latest/samples_report/LUAD_FFPE_Cases.html") %>% 
        html_nodes("table") %>% 
        html_table()

ffpe <- ffpe[[1]] %>% 
        as_tibble()


# Write -------------------------------------------------------------------

write_tsv(ffpe, "./data/tcga-luad/ffpe-cases.tsv")
