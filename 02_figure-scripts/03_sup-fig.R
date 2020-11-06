
# Description -------------------------------------------------------------

# Supplemental Figure 3


# Prepare Workspace --------------------------------------------------

library(tidyverse)
library(DESeq2)
library(survival)
library(survminer)
library(viridis)


# Read in Data ------------------------------------------------------------

imvigor <- read_rds("./data/imvigor210/dds-gsva.Rds")


# Extract colData ---------------------------------------------------------

col_data <- colData(imvigor) %>% 
        as_tibble() %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
        mutate(Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, 
                                                        levels = c("CR", "PR", "SD", "PD", "NE"))) %>% 
        mutate(t_bin = if_else(cd8_t_eff > 0, "hi", "lo"),
               b_bin = if_else(b_cell > 0, "hi", "lo")) %>% 
        unite(b8t, b_bin, t_bin, remove = F) %>% 
        mutate(b8t = factor(b8t, levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo")))


# Fig S3a: PD-L1 IC vs Response --------------------------------------

fig_s3a <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE",
               !is.na(IC.Level)) %>% 
        mutate(IC.Level = factor(IC.Level, levels = c("IC0", "IC1", "IC2+")),
               Best.Confirmed.Overall.Response = factor(Best.Confirmed.Overall.Response, levels = c("CR", "PR", "SD", "PD"))) 

ggplot(fig_s3a, aes(x = IC.Level,
                   fill = Best.Confirmed.Overall.Response)) +
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) +
        theme_minimal() + 
        theme(legend.position = "top",
              text = element_text(size = 15),
              panel.grid = element_blank()) + 
        labs(fill = "", x = NULL, y = "Proportion")
ggsave("./figures/fig_s3/fig_s3a.png", width = 5, height = 7)



# Fig S3b: PD-L1 IC vs B8T ------------------------------------------------

fig_s3b <- col_data %>% 
        filter(!is.na(IC.Level)) %>% 
        mutate(b8t = factor(b8t, 
                            levels = c("hi_hi", "lo_hi", "lo_lo", "hi_lo"),
                            labels = c("Hi/Hi", "Lo/Hi", "Lo/Lo", "Hi/Lo")),
               IC.Level = factor(IC.Level, levels = c("IC2+", "IC1", "IC0")))

ggplot(fig_s3b, aes(x = b8t, fill = IC.Level)) + 
        geom_bar(position = "fill") + 
        scale_fill_viridis_d(end = 0.95, direction = -1) +
        theme_minimal() + 
        labs(y = "Proportion", 
             fill = "IC Level") +
        theme(text = element_text(size = 15),
              legend.position = "top",
              panel.grid = element_blank())
ggsave("./figures/fig_s3/fig_s3b.png", width = 6, height = 7)
