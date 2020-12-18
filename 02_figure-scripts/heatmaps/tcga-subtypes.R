
# Description -------------------------------------------------------------

# Signatures across TCGA Subtypes


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(DESeq2)
library(viridisLite)


# Read in Data ------------------------------------------------------------

blca <- read_rds("./data/tcga-blca/C_04_merge-gsva.rds")

blca_exp <- assay(blca, 2) %>% 
        t() %>% 
        scale() %>% 
        t() %>% 
        apply(2, function(x) ifelse(x < -3, -3, x)) %>% 
        apply(2, function(x) ifelse(x >  3,  3, x))

assay(blca, 2, withDimnames = F) <- blca_exp

colnames(blca) <- blca$sample

signatures <- read_rds("./data/signatures/signatures.Rds") %>% 
        unlist() %>% 
        as_tibble(rownames = "cell") %>% 
        mutate(cell = str_remove(cell, "[:digit:]{1,2}$"))

col_data <- colData(blca) %>% 
        as_tibble() %>% 
        mutate(b_bin = if_else(b_cell > 0, "hi", "lo"))


# Rosenberg CD8 ------------------------------------------------------

col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_rose > 0,  "hi", "lo"))

annotation <- 
        data.frame(row.names = colnames(blca),
                   b_bin = col_data$b_bin, 
                   t_bin = col_data$t_bin,
                   subtype = col_data$mRNA.cluster) %>% 
        arrange(subtype, t_bin, b_bin)

signatures_rose <- signatures %>% 
        filter(cell %in% c("b_cell", "cd8_rose"))

blca_rose <- blca[which(rownames(blca) %in% signatures_rose$value)]

annotation_rows <- 
        data.frame(gene = rownames(blca_rose)) %>% 
        left_join(signatures_rose, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

blca_sig <- blca_rose[,match(rownames(annotation), colnames(blca_rose))]

pheatmap(assay(blca_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows, 
         color = cividis(100), 
         cluster_cols = F, 
         show_colnames = F, gaps_col = c(73, 120, 126, rep(142, 5), 
                                         146, 156, rep(168, 5), 
                                         214, 220, 236, rep(246, 5), 
                                         263, 269, 289, rep(388, 5), 
                                         392, 394))

