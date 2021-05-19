
# Description -------------------------------------------------------------

# ESTIMATE stromal score of imvigor tumors


# Prepare Workspace -------------------------------------------------------

library(estimate)
library(tidyverse)
library(DESeq2)
library(showtext)
library(gt)

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

col_data <- colData(imv) %>%
        as_tibble() %>% 
        dplyr::filter(!is.na(TC.Level)) %>%
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, 
                            levels = c("lo_lo", "lo_hi", "hi_lo", "hi_hi"),
                            labels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi"))) %>%
        dplyr::select(b8t, TC.Level) %>%
        group_by(TC.Level, b8t) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        pivot_wider(names_from = TC.Level, values_from = n)
        summarize()

col_data %>%
        gt(rowname_col = "b8t") %>% 
        gtsave("./figures/fig_s7/fig_s7.png")
