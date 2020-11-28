
# Description -------------------------------------------------------------

# Read in and tidy GSE13507


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(GEOquery)


# Read in Data ------------------------------------------------------------

dat <- getGEO("GSE13507", parseCharacteristics = F, GSEMatrix = T, getGPL = T)
dat <- dat[[1]]


# Tidy Expression Data -----------------------------------------------

expression <- dat@assayData$exprs 
all(rownames(expression) == dat@featureData@data$ID)
rownames(expression) <- dat@featureData@data$ILMN_Gene


# Tidy Clinical Data -------------------------------------------------

clin <- dat@phenoData@data


# Write Files --------------------------------------------------------

write_rds(expression, compress = "gz", file = "./data/kim/01_tidy-expression.Rds")
write_rds(clin, compress = "gz", file = "./data/kim/01_tidy-clin.Rds")
