
# Description -------------------------------------------------------------

# Heatmaps of signature genes throughout tumors


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(DESeq2)
library(viridisLite)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")

colnames(imvigor) <- imvigor$sample

signatures <- read_rds("./data/signatures/signatures.Rds") %>% 
        unlist() %>% 
        as_tibble(rownames = "cell") %>% 
        mutate(cell = str_remove(cell, "[:digit:]{1,2}$"))

imvigor_sig <- imvigor[which(rownames(imvigor) %in% signatures$value)]

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(b_bin = if_else(b_cell > 0, "hi", "lo"),
               t_bin = if_else(cd8_t_eff > 0,  "hi", "lo"))

annotation <- data.frame(row.names = colnames(imvigor), b_bin = col_data$b_bin, t_bin = col_data$t_bin) %>% 
        arrange(b_bin, t_bin)

annotation_rows <- data.frame(gene = rownames(imvigor_sig)) %>% 
        left_join(signatures, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

imvigor_sig <- imvigor_sig[,match(rownames(annotation), colnames(imvigor_sig))]

pheatmap(assay(imvigor_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows,
         scale = "row", color = cividis(100), cluster_cols = F, show_colnames = F)
