
# Description -------------------------------------------------------------

# Figure 3: The B8T high/high signature associates with an overall survival
# advantage post-CPI in men with metastatic bladder cancer post-CPI, but not
# women.


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(viridis)
library(DESeq2)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")


# Extract colData ---------------------------------------------------------

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               gender = factor(gender, levels = c("male", "female")))


# Fig 3a: Sex vs Survival -------------------------------------------------

png(filename = "./figures/fig_3/fig_3a.png", width = 5.5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ gender, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex",
                   legend.labs = c("M", "F"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 3b: Sex vs B8T vs Survival ------------------------------------------

hi_hi <- col_data %>% 
        filter(b8t == "hi_hi")

png(filename = "./figures/fig_3/fig_3b_i.png", width = 5.5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ gender, data = hi_hi) %>% 
        ggsurvplot(pval = T, 
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex\n(B8T = Hi/Hi)",
                   legend.labs = c("Male", "Female"),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


lo_lo <- col_data %>% 
        filter(b8t == "lo_lo")

png(filename = "./figures/fig_3/fig_3b_ii.png", width = 5.5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ gender, data = lo_lo) %>% 
        ggsurvplot(pval = T, 
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex\n(B8T = Lo/Lo)",
                   legend.labs = c("Male", "Female"),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 3c: B8T vs Survival - Males -----------------------------------------

male <- col_data %>% 
        filter(gender == "male")

png(filename = "./figures/fig_3/fig_3c.png", width = 5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ b8t, data = male) %>% 
        ggsurvplot(pval = T, 
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.title = "B8T\n(Sex = Males)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 3d: B8T vs Survival - Females -----------------------------------------

female <- col_data %>% 
        filter(gender == "female", b8t != "lo_hi")

png(filename = "./figures/fig_3/fig_3d.png", width = 5, height = 4, units = "in", res = 288)
survfit(Surv(os, censOS) ~ b8t, data = female) %>% 
        ggsurvplot(pval = T, 
                   palette = c("#DEE318FF", "#32648EFF", "#440154FF"), 
                   legend.title = "B8T\n(Sex = Females)",
                   legend.labs = c("Hi/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 10,
                   legend = "right",
                   pval.coord = c(0, 0.05))
dev.off()


# Fig 3e: Sex vs B8T vs Response ------------------------------------------

labels <- c("Female", "Male")
names(labels) <- c("female", "male")

fig_3e <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

ggplot(fig_3e, aes(b8t, fill = Best.Confirmed.Overall.Response)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) + 
        facet_wrap(~gender,  labeller = labeller(gender = labels)) +
        labs(fill = "", y = "Proportion", x = "B8T") +
        theme_minimal() +
        scale_x_discrete(labels = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo")) +
        theme(text = element_text(size = 15),
              panel.grid = element_blank(),
              strip.text = element_text(size = 15),
              panel.spacing = unit(2, "lines"),
              legend.position = "top")
ggsave("./figures/fig_3/fig_3e.png", width = 6, height = 4)
