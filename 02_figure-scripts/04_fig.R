
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


# Fig 4a: Sex vs Survival -------------------------------------------------

png(filename = "./figures/fig_4/fig_4a.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ gender, data = col_data) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.22,
                   title = "A.",
                   risk.table.title = "No. at risk",
                   palette = c("#1D557D", "#EB9BC7"),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "Sex",
                   legend.labs = c("M", "F"), 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()



# Fig 4b: B8T vs Survival - Males -----------------------------------------

male <- col_data %>% 
        filter(gender == "male")

png(filename = "./figures/fig_4/fig_4b.png", width = 5.5, height = 5.7, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = male) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.3,
                   title = "B.",
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1), 
                   legend.title = "B8T\n(Sex = M)",
                   legend.labs = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 4c: B8T vs Survival - Females -----------------------------------------

female <- col_data %>% 
        filter(gender == "female", b8t != "lo_hi")

png(filename = "./figures/fig_4/fig_4c.png", width = 5.5, height = 5.6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = female) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   title = "C.",
                   risk.table.title = "No. at risk",
                   palette = c("#DEE318FF", "#32648EFF", "#440154FF"), 
                   legend.title = "B8T\n(Sex = F)",
                   legend.labs = c("Hi/Hi", "Lo/Lo", "Hi/Lo"), 
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 4d: TMB vs Sex vs Survival ------------------------------------------

fig_4d <- col_data %>% 
        filter(!is.na(tmb_bins)) %>% 
        unite(s_tmb, gender, tmb_bins, remove = F) %>% 
        mutate(s_tmb = factor(s_tmb, levels=c("male_Hi", "male_Lo", "female_Hi", "female_Lo")))

png(filename = "./figures/fig_4/fig_4d.png", width = 5.5, height = 5.7, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ s_tmb, data = fig_4d) %>% 
        ggsurvplot(pval = T,
                   risk.table = T,
                   risk.table.height = 0.3,
                   title = "D.",
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   legend.title = "Sex/TMB",
                   legend.labs = c("M/Hi", "M/Lo", "F/Hi", "F/Lo"),
                   xlab = "Overall Survival (Months)", 
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 4e: PD-L1 IC vs Survival (Males) ------------------------------------

fig_4e <- col_data %>% 
        filter(gender == "male")

png(filename = "./figures/fig_4/fig_4e.png", width = 5.5, height = 5.6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level, data = fig_4e) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   title = "E.",
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "PD-L1\nIC Level\n(Sex = M)",
                   legend.labs = c("IC2+", "IC1", "IC0"),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()


# Fig 4f: PD-L1 IC vs Survival (Females) ----------------------------------

fig_4f <- col_data %>% 
        filter(gender == "female")

png(filename = "./figures/fig_4/fig_4f.png", width = 5.5, height = 5.6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level, data = fig_4f) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   title = "F.",
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   xlab = "Overall Survival (Months)", 
                   legend.title = "PD-L1\nIC Level\n(Sex = F)",
                   legend.labs = c("IC2+", "IC1", "IC0"),
                   font.xtickslab = 15, 
                   font.ytickslab = 15, 
                   font.legend = 15,
                   fontsize = 5,
                   legend = "right",
                   pval.coord = c(0, 0.05))
ggsurv$plot <- ggsurv$plot + 
        theme(plot.title.position = "plot",
              plot.title = element_text(face = "bold", size = 27),
              axis.title.y = element_text(size = 20),
              axis.title.x = element_text(size = 20))
ggsurv$table <- ggsurv$table +
        ylab(NULL) + 
        xlab(NULL) +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_text(size = 15),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              plot.title = element_text(size = 20))
ggsurv
dev.off()
