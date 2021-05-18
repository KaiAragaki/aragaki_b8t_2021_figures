
# Description -------------------------------------------------------------

# ESTIMATE stromal score of imvigor tumors


# Prepare Workspace -------------------------------------------------------

library(estimate)
library(tidyverse)
library(DESeq2)
library(showtext)

showtext_auto()

if(Sys.info()[["sysname"]] == "Windows") {
        font_add("GillSans", "GIL_____.TTF")
} else if(Sys.info()[["sysname"]] == "Darwin") {
        font_add("GillSans", "GillSans.ttc")
}

theme_tufte <- function(font_size = 30) {
        theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "#FFFFF8", color = "#CCCCCC"),
                plot.background = element_rect(fill = "#FFFFF8"),
                strip.background = element_rect(fill = "#BBBBB0"),
                legend.background = element_rect(fill = "#FFFFF8"), 
                legend.position = "top",
                legend.key = element_blank(),
                text = element_text(family = "GillSans", size = font_size)
        )
}


# Read in Data ------------------------------------------------------------

imv <- read_rds("./data/imvigor210/dds-gsva.Rds")

imv_2 <- imv %>% 
        assay(2) %>%
        as.data.frame()

colnames(imv_2) <- imv$sample


outputGCT(imv_2, "./data/imvigor210/imv.gct")
estimateScore("./data/imvigor210/imv.gct", "./data/imvigor210/imv_est-out.tsv", platform = "illumina")

a <- read_tsv("./data/imvigor210/imv_est-out.tsv")
colnames(a) <- a[2,]

a <- dplyr::slice(a, -(1:2))

a <- t(a)

colnames(a) <- a[1,]
a <- a %>% 
        as.data.frame() %>% 
        dplyr::slice(-(1:2)) %>% 
        mutate(across(everything(), as.numeric))

colData(imv) <- DataFrame(cbind(a, colData(imv)))

col_data <- colData(imv) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               tmb_bins = if_else(FMOne.mutation.burden.per.MB < 10, "Lo", "Hi"),
               tmb_bins = factor(tmb_bins, levels = c("Hi", "Lo"))) %>% 
        unite(tmb_b8t, tmb_bins, b_bin, t_bin, remove = F) %>% 
        mutate(tmb_b8t = factor(tmb_b8t, levels = c("Hi_hi_hi", "Hi_lo_hi", "Hi_lo_lo", "Hi_hi_lo", "Lo_hi_hi", "Lo_lo_hi", "Lo_lo_lo", "Lo_hi_lo")))


ggplot(col_data, aes(x = b8t, y = ImmuneScore, fill = b8t)) +
        scale_color_viridis_d(option = "plasma", end = 0.8) +
        scale_fill_viridis_d(option = "plasma", end = 0.8) + 
        geom_boxplot(alpha = 0.2) + 
        geom_jitter(aes(color = b8t, alpha = 0.8)) + 
        theme_tufte() + 
        theme(legend.position = "none")

ggplot(col_data, aes(x = ImmuneScore, y = StromalScore, color = b8t)) +
        scale_color_viridis_d(option = "plasma", end = 0.8) +
        scale_fill_viridis_d(option = "plasma", end = 0.8) + 
        geom_point() + 
        theme_tufte()
