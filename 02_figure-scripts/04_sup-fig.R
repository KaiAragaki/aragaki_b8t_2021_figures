
# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(viridis)
library(DESeq2)
library(broom)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")


# Extract colData ---------------------------------------------------------

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               tmb_bins = if_else(FMOne.mutation.burden.per.MB < 10, "Lo", "Hi"),
               tmb_bins = factor(tmb_bins, levels = c("Lo", "Hi")),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               gender = factor(gender, levels = c("male", "female")),
               IC.Level = factor(IC.Level, levels = c("IC2+", "IC1", "IC0")))


# Fig S4a: B-Cell or T-cell Signature Distribution across Sex -------------

labels <- c("Male", "Female")
names(labels) <- c("male", "female")

ggplot(col_data, aes(b_cell, fill = gender, color = gender)) +
        geom_density(color = NA) +
        scale_fill_discrete(type = c("#1D557D", "#EB9BC7")) +
        xlab("B-cell Signature") +
        ylab("Density") + 
        labs(fill = "Sex") +
        facet_wrap(~gender, labeller = labeller(gender = labels)) +
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA)) +
        coord_cartesian(xlim = c(-1, 1))
ggsave("./figures/fig_s4/fig_s4a_i.png", width = 6, height = 3)

ggplot(col_data, aes(cd8_rose, fill = gender, color = gender)) +
        geom_density(color = NA) +
        scale_fill_discrete(type = c("#1D557D", "#EB9BC7")) +
        xlab("CD8+ T-cell Signature") +
        ylab("Density") + 
        labs(fill = "Sex") +
        facet_wrap(~gender, labeller = labeller(gender = labels)) +
        theme_minimal() +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA)) +
        coord_cartesian(xlim = c(-1, 1))
ggsave("./figures/fig_s4/fig_s4a_ii.png", width = 6, height = 3)




# Fig s4b: Sex vs B8T vs Response -----------------------------------------

labels <- c("Female", "Male")
names(labels) <- c("female", "male")

fig_s4 <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

ggplot(fig_s4, aes(b8t, fill = Best.Confirmed.Overall.Response)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) + 
        facet_wrap(~gender,  labeller = labeller(gender = labels)) +
        labs(fill = "", y = "Proportion", x = "B-cell/CD8+ T-cell Signature") +
        theme_minimal() +
        scale_x_discrete(labels = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo")) +
        theme(text = element_text(size = 15),
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              panel.spacing = unit(2, "lines"),
              legend.position = "top")
ggsave("./figures/fig_s4/fig_s4b.png", width = 6, height = 4)


# Fig s4c: Sex vs BCGS -----------------------------------------------

male <- filter(col_data, gender == "male")

png(filename = "./figures/fig_s4/fig_s4ci.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b_bin, data = male) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   xlab = "Overall Survival (Months)", palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "BCGS\n(Sex = Male)",
                   legend.labs = c("Hi", "Lo"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

female <- filter(col_data, gender == "female")

png(filename = "./figures/fig_s4/fig_s4cii.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b_bin, data = female) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   xlab = "Overall Survival (Months)", palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "BCGS\n(Sex = Female)",
                   legend.labs = c("Hi", "Lo"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


# Fig s4d: Sex vs CD8GS ----------------------------------------------


male <- filter(col_data, gender == "male")

png(filename = "./figures/fig_s4/fig_s4di.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ t_bin, data = male) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   xlab = "Overall Survival (Months)", palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "CD8TGS\n(Sex = Male)",
                   legend.labs = c("Hi", "Lo"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()

female <- filter(col_data, gender == "female")

png(filename = "./figures/fig_s4/fig_s4dii.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ t_bin, data = female) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   xlab = "Overall Survival (Months)", palette = viridis(2, end = 0.95, direction = -1), 
                   legend.title = "CD8TGS\n(Sex = Female)",
                   legend.labs = c("Hi", "Lo"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank())
ggsurv
dev.off()


