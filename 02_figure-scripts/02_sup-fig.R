
# Description -------------------------------------------------------------

# Supplemental Figure 2: Patients with TMB (high) tumors have higher 
# proportion of responders


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
        mutate(t_bin = if_else(cd8_rose > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo"),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")),
               tmb_bins = if_else(FMOne.mutation.burden.per.MB < 10, "Lo", "Hi"),
               tmb_bins = factor(tmb_bins, levels = c("Lo", "Hi"))) %>% 
        unite(tmb_b8t, tmb_bins, b_bin, t_bin, remove = F) %>% 
        mutate(tmb_b8t = factor(tmb_b8t, levels = c("Hi_hi_hi", "Hi_lo_hi", "Hi_lo_lo", "Hi_hi_lo", "Lo_hi_hi", "Lo_lo_hi", "Lo_lo_lo", "Lo_hi_lo")))

# Supplemental Fig 2a --------------------------------------------------

fig_s2a <- filter(col_data, Best.Confirmed.Overall.Response != "NE", !is.na(tmb_bins))

ggplot(fig_s2a, aes(tmb_bins, fill = Best.Confirmed.Overall.Response)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) + 
        theme_minimal() + 
        theme(legend.position = "top",
              text = element_text(size = 15),
              panel.grid = element_blank()) + 
        labs(fill = "", x = NULL, y = "Proportion") +
        scale_x_discrete(labels = c("TMB Lo", "TMB Hi"))
ggsave("./figures/fig_s2/fig_s2a.png", width = 4, height = 3)


# Supplemental Fig 2b --------------------------------------------------

fig_s2b <- col_data %>% 
        filter(!is.na(tmb_bins)) 

labels <- c("TMB Lo", "TMB Hi")
names(labels) <- c("Lo", "Hi")

table_fig_s2b <- fig_s2b %>%
        group_by(tmb_bins, t_bin, b_bin) %>% 
        summarize(num = n()) %>%
        group_by(tmb_bins) %>% 
        mutate(prop = num/sum(num))

ggplot(table_fig_s2b, aes(x = b_bin, y = t_bin, color = as.factor(rep(1:2, each = 4)))) + 
        geom_point(aes(size = prop)) + 
        scale_color_viridis_d(begin = 0.1, end = 0.5) + 
        facet_wrap(~tmb_bins, labeller = labeller(tmb_bins = labels)) + 
        theme_minimal() + 
        scale_size_area(max_size = 45) + 
        geom_text(aes(label = round(prop, 2)), color = "white",  size = 7) +
        theme(text = element_text(size = 15),
              legend.position = "none",
              panel.grid = element_blank(),
              panel.spacing = unit(2, "lines"),
              panel.background = element_rect(fill = "gray95", color = NA), 
              strip.text = element_text(size = 15)) +
        labs(x = "B-cell Signature", y = "CD8+ T-cell Signature") +
        scale_x_discrete(labels = c("Lo", "Hi")) +
        scale_y_discrete(labels = c("Lo", "Hi"))
ggsave("./figures/fig_s2/fig_s2b.png", width = 6, height = 3.5)

# Supplemental Fig 2c --------------------------------------------------

fig_s2c <- col_data %>% 
        filter(!is.na(tmb_bins))

mt <- pairwise_survdiff(Surv(os, censOS) ~ tmb_b8t, data = fig_s2c)$p.value

fig_s2c <- fig_s2c %>% 
        filter(b8t == "lo_lo")

png(filename = "./figures/fig_s2/fig_s2c.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ tmb_b8t, data = fig_s2c) %>% 
        ggsurvplot(pval = round(mt["Lo_lo_lo", "Hi_lo_lo"], digits = 4), 
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   legend.title = "TMB\n(B8T = Lo/Lo)",
                   legend.labs = c("Hi", "Lo"),
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
