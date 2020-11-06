
# Description -------------------------------------------------------------

# Figure 4: The B8T high/high signature, or high PD-L1 immune cell levels,
# associates with higher OS in response to anti-PD-L1 therapy in men with
# metastatic urothelial cancer, but not women.


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
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               tmb_bins = if_else(FMOne.mutation.burden.per.MB < 10, "Lo", "Hi"),
               tmb_bins = factor(tmb_bins, levels = c("Lo", "Hi")),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               gender = factor(gender, levels = c("male", "female")),
               IC.Level = factor(IC.Level, levels = c("IC2+", "IC1", "IC0")))


# Fig 4a: Sex vs Survival -------------------------------------------------

png(filename = "./figures/fig_4/fig_4a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ gender, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex",
                   legend.labs = c("M", "F"), 
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


# Fig 4b: Sex vs B8T vs Survival ------------------------------------------

hi_hi <- col_data %>% 
        filter(b8t == "hi_hi")

mt <- pairwise_survdiff(Surv(os, censOS) ~ b8t + gender, data = col_data, p.adjust.method = "none") %>% 
        tidy()

considered_comparisons <- mt %>% 
        separate(group1, c("b8t_1", "sex_1"), sep = ",") %>% 
        separate(group2, c("b8t_2", "sex_2"), sep = ",") %>% 
        filter(b8t_1 == b8t_2) %>% 
        mutate(adj = p.adjust(p.value, "BH"))

png(filename = "./figures/fig_4/fig_4b.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ gender, data = hi_hi) %>% 
        ggsurvplot(pval = round(considered_comparisons$adj[1], 2),
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex\n(B8T = Hi/Hi)",
                   legend.labs = c("Male", "Female"),
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


# Fig 4c: Sex vs B8T vs Survival ------------------------------------------

lo_lo <- col_data %>% 
        filter(b8t == "lo_lo")

mt <- pairwise_survdiff(Surv(os, censOS) ~ b8t + gender, data = col_data, p.adjust.method = "none") %>% 
        tidy()

considered_comparisons <- mt %>% 
        separate(group1, c("b8t_1", "sex_1"), sep = ",") %>% 
        separate(group2, c("b8t_2", "sex_2"), sep = ",") %>% 
        filter(b8t_1 == b8t_2) %>% 
        mutate(adj = p.adjust(p.value, "BH"))

png(filename = "./figures/fig_4/fig_4c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ gender, data = lo_lo) %>% 
        ggsurvplot(pval = round(considered_comparisons$adj[3], 2),
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex\n(B8T = Lo/Lo)",
                   legend.labs = c("Male", "Female"),
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



# Fig 4d: B8T vs Survival - Males -----------------------------------------

male <- col_data %>% 
        filter(gender == "male")

png(filename = "./figures/fig_4/fig_4d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = male) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.title = "B8T\n(Sex = Males)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
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


# Fig 4e: B8T vs Survival - Females -----------------------------------------

female <- col_data %>% 
        filter(gender == "female", b8t != "lo_hi")

png(filename = "./figures/fig_4/fig_4e.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = female) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   risk.table.title = "No. at risk",
                   palette = c("#DEE318FF", "#32648EFF", "#440154FF"), 
                   legend.title = "B8T\n(Sex = Females)",
                   legend.labs = c("Hi/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
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


# Fig 4f: TMB vs Sex vs Survival ------------------------------------------

fig_4f <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        unite(s_tmb, gender, tmb_bins) %>% 
        mutate(s_tmb = factor(s_tmb, levels=c("male_Hi", "male_Lo", "female_Hi", "female_Lo")))

png(filename = "./figures/fig_4/fig_4f.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ s_tmb, data = fig_4f) %>% 
        ggsurvplot(pval = T,
                   risk.table = T,
                   risk.table.height = 0.3,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "Sex/TMB",
                   legend.labs = c("M/Hi", "M/Lo", "F/Hi", "F/Lo"),
                   xlab = "Overall Survival (Months)", 
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



# Fig 4g: PD-L1 IC vs Survival (Males) ------------------------------------

fig_4g <- col_data %>% 
        filter(gender == "male")

png(filename = "./figures/fig_4/fig_4g.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level, data = fig_4g) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "IC Level\n(Sex = M)",
                   legend.labs = c("IC2+", "IC1", "IC0"),
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


# Fig 4h: PD-L1 IC vs Survival (Females) ----------------------------------

fig_4h <- col_data %>% 
        filter(gender == "female")

png(filename = "./figures/fig_4/fig_4h.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level, data = fig_4h) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "IC Level\n(Sex = F)",
                   legend.labs = c("IC2+", "IC1", "IC0"),
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

