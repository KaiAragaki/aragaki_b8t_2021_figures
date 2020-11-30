
# Description -------------------------------------------------------------

# Read in and tidy GSE48075


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(GEOquery)


# Read in Data ------------------------------------------------------------

dat <- getGEO("GSE48075", parseCharacteristics = F, GSEMatrix = T, getGPL = T)
dat <- dat[[1]]


# Tidy Annotation Data -----------------------------------------------

annot <- data.frame(symbol = dat@featureData@data$Symbol, 
                    nuID = dat@featureData@data$nuID)


# Tidy Expression Data -----------------------------------------------

expression <- dat@assayData$exprs 
all(rownames(expression) == dat@featureData@data$ID)
rownames(expression) <- dat@featureData@data$nuID
expression <- log2(expression)
expression <- sweep(expression, 1, apply(expression, 1, median))


# Tidy Clinical Data -------------------------------------------------

clin <- dat@phenoData@data

# Sex data received from Woonyoung Choi 2020-11-24
sex <- read_tsv("./data/choi/142_gender.txt")

clin <- inner_join(clin, sex, by = c("title" = "X1"))


# Write Files --------------------------------------------------------

write_rds(expression, compress = "gz", file = "./data/choi/01_tidy-expression.Rds")
write_rds(clin, compress = "gz", file = "./data/choi/01_tidy-clin.Rds")
write_rds(annot, compress = "gz", file = "./data/choi/01_tidy-annot.Rds")