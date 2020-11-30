
# Description -------------------------------------------------------------

# Read in and tidy GSE32894


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(GEOquery)


# Read in Data ------------------------------------------------------------

dat <- getGEO("GSE32894", parseCharacteristics = F, GSEMatrix = T, getGPL = T)
dat <- dat[[1]]


# Tidy Annotation Data -----------------------------------------------

annot <- data.frame(symbol = dat@featureData@data$Symbol, 
                    nuID = dat@featureData@data$nuID)


# Tidy Expression Data -----------------------------------------------

expression <- dat@assayData$exprs 
all(rownames(expression) == dat@featureData@data$ID)
rownames(expression) <- dat@featureData@data$nuID


# Tidy Clinical Data -------------------------------------------------

clin <- dat@phenoData@data


# Write Files --------------------------------------------------------

write_rds(expression, compress = "gz", file = "./data/sjodahl/01_tidy-expression.Rds")
write_rds(clin, compress = "gz", file = "./data/sjodahl/01_tidy-clin.Rds")
write_rds(annot, compress = "gz", file = "./data/sjodahl/01_tidy-annot.Rds")
