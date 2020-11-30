
# Description -------------------------------------------------------------

# Read in and tidy GSE13507


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(GEOquery)
library(lumi)


# Read in Data ------------------------------------------------------------

dat <- getGEO("GSE13507", parseCharacteristics = F, GSEMatrix = T, getGPL = T)
dat <- dat[[1]]


# Tidy Annotation Data -----------------------------------------------

annot <- data.frame(symbol = dat@featureData@data$Symbol, 
                    nuID = seq2id(dat@featureData@data$SEQUENCE))


# Tidy Expression Data -----------------------------------------------

expression <- dat@assayData$exprs 
all(rownames(expression) == dat@featureData@data$ID)
rownames(expression) <- annot$nuID
expression <- sweep(expression, 1, apply(expression, 1, median))


# Tidy Clinical Data -------------------------------------------------

clin <- dat@phenoData@data


# Write Files --------------------------------------------------------

write_rds(expression, compress = "gz", file = "./data/kim/01_tidy-expression.Rds")
write_rds(clin, compress = "gz", file = "./data/kim/01_tidy-clin.Rds")
write_rds(annot, compress = "gz", file = "./data/kim/01_tidy-annot.Rds")