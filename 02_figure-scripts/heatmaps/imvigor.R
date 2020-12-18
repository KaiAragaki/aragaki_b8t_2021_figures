
# Description -------------------------------------------------------------

# Heatmaps of signature genes throughout tumors


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(DESeq2)
library(viridisLite)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")

# Compress Z score to preserve dynamic range

imvigor_exp <- assay(imvigor, 2) %>% 
        t() %>% 
        scale() %>% 
        t() %>% 
        apply(2, function(x) ifelse(x < -3, -3, x)) %>% 
        apply(2, function(x) ifelse(x >  3,  3, x))

assay(imvigor, 2, withDimnames = F) <- imvigor_exp

colnames(imvigor) <- imvigor$sample

signatures <- read_rds("./data/signatures/signatures.Rds") %>% 
        unlist() %>% 
        as_tibble(rownames = "cell") %>% 
        mutate(cell = str_remove(cell, "[:digit:]{1,2}$"))

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(b_bin = if_else(b_cell > 0, "hi", "lo"))


# Rosenberg CD8 ------------------------------------------------------

col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_rose > 0,  "hi", "lo"))

annotation <- 
        data.frame(row.names = colnames(imvigor), 
                   b_bin = col_data$b_bin, 
                   t_bin = col_data$t_bin) %>% 
        arrange(b_bin, t_bin)

signatures_rose <- signatures %>% 
        filter(cell %in% c("b_cell", "cd8_rose"))

imvigor_rose <- imvigor[which(rownames(imvigor) %in% signatures_rose$value)]

annotation_rows <- 
        data.frame(gene = rownames(imvigor_rose)) %>% 
        left_join(signatures_rose, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

imvigor_sig <- imvigor_rose[,match(rownames(annotation), colnames(imvigor_rose))]

pheatmap(assay(imvigor_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows, color = cividis(100), cluster_cols = F, show_colnames = F)


# Prat CD8 ------------------------------------------------------

col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_prat > 0,  "hi", "lo"))

annotation <- 
        data.frame(row.names = colnames(imvigor), 
                   b_bin = col_data$b_bin, 
                   t_bin = col_data$t_bin) %>% 
        arrange(b_bin, t_bin)

signatures_prat <- signatures %>% 
        filter(cell %in% c("b_cell", "cd8_prat"))

imvigor_prat <- imvigor[which(rownames(imvigor) %in% signatures_prat$value)]

annotation_rows <- 
        data.frame(gene = rownames(imvigor_prat)) %>% 
        left_join(signatures_prat, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

imvigor_sig <- imvigor_prat[,match(rownames(annotation), colnames(imvigor_prat))]

pheatmap(assay(imvigor_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows, color = cividis(100), cluster_cols = F, show_colnames = F)


# Fehr CD8 ------------------------------------------------------

col_data <- col_data %>% 
        mutate(t_bin = if_else(cd8_fehr > 0,  "hi", "lo"))

annotation <- 
        data.frame(row.names = colnames(imvigor), 
                   b_bin = col_data$b_bin, 
                   t_bin = col_data$t_bin) %>% 
        arrange(b_bin, t_bin)

signatures_fehr <- signatures %>% 
        filter(cell %in% c("b_cell", "cd8_fehr"))

imvigor_fehr <- imvigor[which(rownames(imvigor) %in% signatures_fehr$value)]

annotation_rows <- 
        data.frame(gene = rownames(imvigor_fehr)) %>% 
        left_join(signatures_fehr, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

imvigor_sig <- imvigor_fehr[,match(rownames(annotation), colnames(imvigor_fehr))]

pheatmap(assay(imvigor_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows, color = cividis(100), cluster_cols = F, show_colnames = F)


# IFNg ---------------------------------------------------------------

col_data <- col_data %>% 
        mutate(ifng_bin = if_else(ifng > 0,  "hi", "lo"))

annotation <- 
        data.frame(row.names = colnames(imvigor), 
                   b_bin = col_data$b_bin, 
                   g_bin = col_data$ifng_bin) %>% 
        arrange(b_bin, g_bin)

signatures_ifng <- signatures %>% 
        filter(cell %in% c("b_cell", "ifng"))

imvigor_ifng <- imvigor[which(rownames(imvigor) %in% signatures_ifng$value)]

annotation_rows <- 
        data.frame(gene = rownames(imvigor_ifng)) %>% 
        left_join(signatures_ifng, by = c("gene" = "value")) %>% 
        column_to_rownames("gene")

imvigor_sig <- imvigor_ifng[,match(rownames(annotation), colnames(imvigor_ifng))]

pheatmap(assay(imvigor_sig, 2), 
         annotation_col = annotation, 
         annotation_row = annotation_rows, color = cividis(100), cluster_cols = F, show_colnames = F)