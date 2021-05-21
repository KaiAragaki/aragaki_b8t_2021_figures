
# Description -------------------------------------------------------------

# ESTIMATE stromal score of IMvigor210 and TCGA BLCA tumors


# Prepare Workspace -------------------------------------------------------

library(estimate)
library(tidyverse)
library(DESeq2)
library(showtext)
library(broom)

showtext_auto()

if(Sys.info()[["sysname"]] == "Windows") {
        font_add("GillSans", "GIL_____.TTF")
} else if(Sys.info()[["sysname"]] == "Darwin") {
        font_add("GillSans", "GillSans.ttc")
}

theme_tufte <- function(font_size = 30) {
        theme(
                panel.grid = element_blank(),
                panel.background = element_rect(fill = "#FFFFF8", 
                                                color = "#CCCCCC"),
                plot.background = element_rect(fill = "#FFFFF8"),
                strip.background = element_rect(fill = "#BBBBB0"),
                legend.background = element_rect(fill = "#FFFFF8"), 
                legend.position = "top",
                legend.key = element_blank(),
                text = element_text(family = "GillSans", 
                                    size = font_size)
        )
}


# Read in Data -----------------------------------------------------------------

imv <- read_rds("./data/imvigor210/dds-gsva.Rds")
blca <- read_rds("./data/tcga-blca/C_04_merge-gsva.rds")


# Convert to GCT ---------------------------------------------------------------

imv_2 <- imv %>% 
        assay(2) %>%
        as.data.frame() %>%
        setNames(imv$sample)

outputGCT(imv_2, "./data/imvigor210/imv.gct")

blca_2 <- blca %>%
        assay(2) %>%
        as.data.frame() %>%
        setNames(blca$sample)


outputGCT(blca_2, "./data/tcga-blca/blca.gct")


# ESTIMATE ---------------------------------------------------------------------

estimateScore("./data/imvigor210/imv.gct", 
              "./data/imvigor210/imv_est-out.tsv", 
              platform = "illumina")

estimateScore("./data/tcga-blca/blca.gct",
              "./data/tcga-blca/blca_est-out.tsv",
              platform = "illumina")

format_coldata <- function(est_path, dds) {
        df <- read_tsv(est_path)
        colnames(df) <- df[2,]
        df <- dplyr::slice(df, -(1:2)) %>%
                t()
        colnames(df) <- df[1,]
        df <- df %>%
                as.data.frame() %>%
                dplyr::slice(-(1:2)) %>%
                mutate(across(everything(), as.numeric))
        print(all(rownames(df) == rownames(colData(dds))))
        df <- df %>%        
                cbind(colData(dds)) %>%
                as_tibble() %>% 
                mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
                       b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
                unite(b8t, b_bin, t_bin, remove = F) %>% 
                mutate(b8t = factor(b8t, 
                                    levels = c("lo_lo", "lo_hi", "hi_lo", "hi_hi"),
                                    labels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi")))
}

cd_imv <- format_coldata("./data/imvigor210/imv_est-out.tsv", imv)
cd_blca <- format_coldata("./data/tcga-blca/blca_est-out.tsv", blca)

# AOV --------------------------------------------------------------------------

aov(StromalScore ~ b8t, data = cd_imv) %>%
        summary()

aov(StromalScore ~ b8t, data = cd_blca) %>%
        summary()


# T-tests -----------------------------------------------------------------

make_tt_df <- function(col_data) {
        col_data %>%
                dplyr::select(StromalScore, b8t) %>%
                pivot_wider(names_from = b8t, values_from = StromalScore, values_fn = list)
}

test_groups <- function(group_1, group_2, tt_df){
        t.test(unlist(tt_df[[group_1]]), unlist(tt_df[[group_2]])) %>%
                tidy() %>%
                mutate(group_1 = group_1,
                       group_2 = group_2)
}

imv_tt <- make_tt_df(cd_imv)
blca_tt <- make_tt_df(cd_blca)

do_ttests <- function(tt_df, col_data) {
        ttests_imv <- bind_rows(
                test_groups("Lo/Hi", "Lo/Lo", tt_df),
                test_groups("Lo/Lo", "Hi/Lo", tt_df),
                test_groups("Lo/Lo", "Hi/Hi", tt_df),
                test_groups("Hi/Hi", "Hi/Lo", tt_df),
                test_groups("Hi/Hi", "Lo/Hi", tt_df),
                test_groups("Hi/Lo", "Lo/Hi", tt_df)
        ) %>%
                mutate(p_adj = p.adjust(p.value, method = "BH"),
                       p_star = case_when(p_adj < 0.001 ~ "***",
                                          p_adj < 0.01 ~ "**",
                                          p_adj < 0.05 ~ "*",
                                          T ~ "NS"),
                       group_1 = factor(group_1, levels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi")),
                       group_2 = factor(group_2, levels = c("Lo/Lo", "Lo/Hi", "Hi/Lo", "Hi/Hi")),
                       g1n = as.numeric(group_1),
                       g2n = as.numeric(group_2),
                       x = (g1n + g2n)/2,
                       length = abs(g2n - g1n),
                       y = max(col_data$StromalScore) + ((length-1)*.3*max(col_data$StromalScore)) + ((g1n)*.03*max(col_data$StromalScore)),
                       y_text = if_else(p_star == "NS", y + .06*max(col_data$StromalScore), y + .03*max(col_data$StromalScore)))
}

imv_ttests <- do_ttests(imv_tt, cd_imv)
blca_ttests <- do_ttests(blca_tt, cd_blca)

plot_stromal_score <- function(col_data, ttests, filename, plot_title) {
        ggplot(col_data, aes(x = b8t, y = StromalScore)) +
                scale_color_viridis_d(option = "plasma", end = 0.8) +
                scale_fill_viridis_d(option = "plasma", end = 0.8) + 
                geom_boxplot(alpha = 0.2, aes(fill = b8t)) + 
                geom_jitter(aes(color = b8t), alpha = 0.5, width = 0.2, shape = 16) + 
                geom_segment(data = ttests, aes(x = group_1, xend = group_2, y = y, yend = y), color = "black") +
                geom_text(data = ttests, aes(x = x, y = y_text, label = p_star), size = 12) +
                theme_tufte(40) + 
                coord_cartesian(ylim = c(0, NA)) +
                labs(x = "B8T Signature", title = plot_title) + 
                theme(legend.position = "none")
        ggsave(paste0("./figures/fig_s6/", filename), width = 2.5, height = 3) 
}

plot_stromal_score(cd_imv, imv_ttests, "fig_s6a.png", "IMvigor210")
plot_stromal_score(cd_blca, blca_ttests, "fig_s6b.png", "TCGA BLCA")
