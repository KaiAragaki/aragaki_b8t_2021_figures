
# Description -------------------------------------------------------------

# Read in and tidy Liu et al. data


# Prepare Workspace -------------------------------------------------------

library(readr)
library(dplyr)
library(readxl)
library(tibble)
library(SummarizedExperiment)


# Fetch Data --------------------------------------------------------------

if(!file.exists("./data/liu/tpm.txt")){
        download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0654-5/MediaObjects/41591_2019_654_MOESM3_ESM.txt", 
                      "./data/liu/tpm.txt")
}

if(!file.exists("./data/liu/clin.xlsx")){
        download.file("https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0654-5/MediaObjects/41591_2019_654_MOESM4_ESM.xlsx", 
                      "./data/liu/clin.xlsx", 
                      mode = "wb")
}

# Read in Data and Tidy ---------------------------------------------------

liu <- read_tsv("./data/liu/tpm.txt") %>% 
        t() %>% 
        `colnames<-`(.[1,])
liu <- liu[-1,]

clin <- read_excel("./data/liu/clin.xlsx", skip = 2, n_max = 145)

clin <- clin[-1,]

liu <- liu[,which(colnames(liu) %in% clin$...1)]

liu2 <- apply(liu, 2, as.numeric)

rownames(liu2) <- rownames(liu)

clin <- clin[which(clin$...1 %in% colnames(liu2)),]

clin <- clin[match(colnames(liu2), clin$...1),]

se <- SummarizedExperiment(liu2, colData = DataFrame(clin), rowData = DataFrame(gene = rownames(liu2)))

assay(se, 2) <- apply(assay(se), 2, function(x) log2(x + 1))

assayNames(se) <- c("tmp", "log_2")

# Write SummarizedExperiment ----------------------------------------------

write_rds(se, "./data/liu/se.Rds", compress = "gz")
