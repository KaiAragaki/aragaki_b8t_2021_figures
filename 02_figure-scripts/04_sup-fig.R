
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



# Fig s4: Sex vs B8T vs Response ------------------------------------------

labels <- c("Female", "Male")
names(labels) <- c("female", "male")

fig_s4 <- col_data %>% 
        filter(Best.Confirmed.Overall.Response != "NE")

ggplot(fig_s4, aes(b8t, fill = Best.Confirmed.Overall.Response)) + 
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
ggsave("./figures/fig_s4/fig_s4.png", width = 6, height = 4)
