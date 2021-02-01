
# Description -------------------------------------------------------------

# Figure 3: The B8T signature stratifies patents with different PD-L1 immune
# cell staining into cohorts with different OS in response to anti-PD-L1
# therapy.


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(survival)
library(survminer)
library(viridis)
library(DESeq2)
library(broom)
library(VennDiagram)



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


# Fig 3a: Survival across PD-L1 IC ----------------------------------------

fig_3a <- col_data %>% 
        filter(!is.na(IC.Level))

png(filename = "./figures/fig_3/fig_3a.png", width = 5.5, height = 5.6, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level, data = fig_3a) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.25,
                   risk.table.title = "No. at risk",
                   palette = viridis(3, end = 0.95, direction = -1),
                   title = "A.",
                   xlab = "Overall Survival (Months)", 
                   legend.title = "PD-L1\nIC Level",
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


# Fig 3b: Survival of B8T across PDL1 2+ ----------------------------------

fig_3b <- col_data %>% 
        filter(!is.na(IC.Level), IC.Level == "IC2+")

png(filename = "./figures/fig_3/fig_3b.png", width = 5.7, height = 5.8, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC.Level + b8t, data = fig_3b) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.27,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   title = "B.",
                   legend.title = "B8T Sig.\n(PD-L1 IC = 2+)",
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


# Fig 3c: Survival of B8T across PDL1 0/1 ---------------------------------

fig_3c <- col_data %>% 
        filter(!is.na(IC.Level), IC.Level %in% c("IC0", "IC1"))

png(filename = "./figures/fig_3/fig_3c.png", width = 5.77, height = 5.8, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ b8t, data = fig_3c) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.27,
                   risk.table.title = "No. at risk",
                   palette = viridis(4, end = 0.95, direction = -1),
                   title = "C.",
                   legend.title = "B8T Sig.\n(PL-L1 IC = 0/1)",
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


# Fig 3d: Comparison of B8T Hi/Hi's with IC Level -------------------------

fig_3d <- col_data %>% 
        mutate(IC_bi = case_when(IC.Level == "IC2+" ~ "IC2+",
                                 IC.Level == "IC1" ~ "IC0/1",
                                 IC.Level == "IC0" ~ "IC0/1"),
               IC_bi = factor(IC_bi, levels = c("IC2+", "IC0/1"))) %>% 
        unite(IC_bi_b8t, IC_bi, b8t, sep = ".", remove = F)

fig_3d_his <- filter(fig_3d, b8t == "hi_hi")

mt <- pairwise_survdiff(Surv(os, censOS) ~ IC_bi_b8t, data = fig_3d, p.adjust.method = "none") %>% 
        tidy() %>% 
        separate(group1, c("ic_1", "b8t_1"), sep = "\\.") %>% 
        separate(group2, c("ic_2", "b8t_2"), sep = "\\.") %>% 
        filter(ic_1 != ic_2) %>%
        filter(b8t_1 == b8t_2) %>% 
        filter(ic_1 != "NA")

adj <- p.adjust(mt$p.value, method = "BH")


png(filename = "./figures/fig_3/fig_3d.png", width = 5.5, height = 5.5, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ IC_bi + b8t, data = fig_3d_his) %>% 
        ggsurvplot(pval = round(adj[1], 3),
                   risk.table = T,
                   risk.table.height = 0.22,
                   risk.table.title = "No. at risk",
                   palette = viridis(2, end = 0.95, direction = -1),
                   title = "D.",
                   legend.title = "PL-L1\nIC Level\n(B8T = Hi/Hi)",
                   xlab = "Overall Survival (Months)", 
                   legend.labs = c("IC2+", "IC0/1"), 
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


# Fig 3e: Comparison of all signatures -------------------------------

hi_hi <- col_data %>% 
        filter(b8t == "hi_hi") %>% 
        mutate(type = "B8T Hi/Hi")

tmb_hi <- col_data %>% 
        filter(tmb_bins == "Hi") %>% 
        mutate(type = "TMB Hi")

ic2 <- col_data %>% 
        filter(IC.Level == "IC2+") %>% 
        mutate(type = "PD-L1 IC2+")

hihihi <- col_data %>% 
        filter(tmb_bins == "Hi", b8t == "hi_hi") %>% 
        mutate(type = "B8T, TMB Hi")

tmb_hi_pdl1_hi <- col_data %>% 
        filter(tmb_bins == "Hi", IC.Level == "IC2+") %>% 
        mutate(type = "TMB Hi, PD-L1 IC2+")

hihihi2 <- col_data %>% 
        filter(IC.Level == "IC2+", b8t == "hi_hi") %>% 
        mutate(type = "B8T Hi, PD-L1 IC2+")

hihihihi <- col_data %>% 
        filter(tmb_bins == "Hi", b8t == "hi_hi", IC.Level == "IC2+") %>% 
        mutate(type = "B8T, TMB Hi, PD-L1 IC2+")        

all_type <- rbind(hi_hi, tmb_hi, ic2, tmb_hi_pdl1_hi, hihihi, hihihi2, hihihihi) %>% 
        mutate(type = factor(type, levels = c("B8T, TMB Hi, PD-L1 IC2+", "B8T, TMB Hi", "TMB Hi, PD-L1 IC2+", "B8T Hi, PD-L1 IC2+", "B8T Hi/Hi", "PD-L1 IC2+", "TMB Hi")))

png(filename = "./figures/fig_3/fig_3e.png", width = 5.8, height = 7, units = "in", res = 288)
ggsurv <- survfit(Surv(os, censOS) ~ type, data = all_type) %>% 
        ggsurvplot(pval = T, 
                   risk.table = T,
                   risk.table.height = 0.4, 
                   risk.table.title = "No. at risk",
                   tables.y.text = F,
                   title = "E.",
                   palette = c("#FF0000FF", "#FDE725FF", "#6BCD5AFF", 
                               "#1E9C89FF", "#31668EFF", "#482576FF", 
                               "#000000FF"),
                   legend.title = "Signature", 
                   legend.labs = c("\nB8T Hi\nTMB Hi\nPD-L1 IC2+\n", 
                                   "\nB8T Hi\nTMB Hi\n", 
                                   "\nTMB Hi\nPD-L1 IC2+\n", 
                                   "\nB8T Hi\nPD-L1 IC2+\n", 
                                   "B8T Hi", 
                                   "PD-L1 IC2+", 
                                   "TMB Hi"),
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


# Fig 3f: Venn Diagram Counts ---------------------------------------------

best <- col_data %>% 
        filter(tmb_bins == "Hi"| IC.Level == "IC2+"| b8t == "hi_hi") %>%
        filter(!is.na(tmb_bins)) %>% 
        select(tmb_bins, IC.Level, b8t) 

hihi <- best %>% 
        filter(tmb_bins != "Hi",
               IC.Level != "IC2+")

tmb_hi <- best %>% 
        filter(b8t != "hi_hi",
               IC.Level != "IC2+")

pdl1_ic2 <- best %>% 
        filter(b8t != "hi_hi",
               tmb_bins != "Hi")

pdl1_ic2_tmbhi <- best %>% 
        filter(IC.Level == "IC2+", 
               tmb_bins == "Hi", 
               b8t != "hi_hi")

pdl1_ic2_hihi <- best %>% 
        filter(IC.Level == "IC2+", 
               tmb_bins != "Hi", 
               b8t == "hi_hi")

hihi_tmb_hi <- best %>% 
        filter(b8t == "hi_hi", 
               IC.Level != "IC2+", 
               tmb_bins == "Hi")

hihi_tmb_hi_hi <- best %>% 
        filter(b8t == "hi_hi", 
               IC.Level == "IC2+", 
               tmb_bins == "Hi")

###Venn diagram for b8t-hi/hi, IC2+ and TMB-Hi
fig_3f <- list(tmb_hi=col_data%>%filter(tmb_bins=="Hi")%>%filter(!is.na(tmb_bins))%>%select(sample)%>%unlist,
               IC2=col_data%>%filter(IC.Level=="IC2+")%>%filter(!is.na(tmb_bins))%>%select(sample)%>%unlist,
               b8t_hi_hi=col_data%>%filter(b8t=="hi_hi")%>%filter(!is.na(tmb_bins))%>%select(sample)%>%unlist)

venn.plot <- venn.diagram( fig_3f,
              height = 5.5, width = 5.5, resolution = 300, units='in',
              compression = "lzw",
              filename = "./figures/fig_3/fig_3f.png",
              output=TRUE,
              category.names = c("TMB Hi", "IC2+", "B8T Hi/Hi"),
              col=c("#440154ff", '#21908dff', '#c7e020ff'),
              fill = c(alpha("#440154ff",0.8), alpha('#21908dff',0.8), alpha('#c7e020ff',0.8)),
              main="F.",
              main.col='black',
              main.cex=2,
              main.fontface = "bold",
              main.pos=c(0.1,1),
              lty = "blank",
              lwd=1,
              cex = 2,
              fontface = "bold",
              cat.cex = 2,
              cat.fontface = "bold",
              cat.col = "black",#c("#440154ff", '#21908dff', '#c7e020ff'),
              cat.pos = c(90, -90,0),
              cat.dist = c(-0.08,-0.08,-0.08),
              rotation = 1,
              rotation.degree = 180,
              margin=0)



