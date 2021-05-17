
# Description -------------------------------------------------------------

# Tidy clinical annotation files

# Clinical annotation (.xml files) downloaded from EBI (2020-09-09)


# Prepare Workspace -------------------------------------------------------

library(tidyverse)
library(xml2)
library(sharepointr)

token <- sharepoint_token()


# Read in Data -------------------------------------------------------

imvigor_clin <- sharepoint_get("https://livejohnshopkins.sharepoint.com/sites/GBCIStorage/Shared%20Documents/Datasets/ImVigor210/samples.zip", token)
file.rename(imvigor_clin, "./data/imvigor210/clin.zip")
unzip("./data/imvigor210/clin.zip", exdir = "./data/imvigor210")


# Tidy Data ---------------------------------------------------------------

parse_dat <- function(path) {
        dat <- read_xml(path) %>% 
                xml_children() %>% 
                xml_children()
        dat <- dat[5] %>% 
                as_list() %>% 
                unlist()
        tib <- tibble(items = dat,              
                      tag = rep(c(1, 0), times = length(dat)/2),
                      val = rep(c(0, 1), times = length(dat)/2))
        tag <- filter(tib, as.logical(tag)) %>% select(items) %>% rename(tag = items)
        val <- filter(tib, as.logical(val)) %>% select(items) %>% rename(val = items)
        tibble(tag, val)
}

parse_dat_name <- function(path) {
        dat <- read_xml(path) %>% 
                xml_children() %>% 
                xml_children()
        dat <- dat[1] %>% 
                as_list() %>% 
                unlist()
        tib <- tibble(items = dat,              
                      tag = rep(c(1, 0), times = length(dat)/2),
                      val = rep(c(0, 1), times = length(dat)/2))
        tag <- filter(tib, as.logical(tag)) %>% select(items) %>% rename(tag = items)
        val <- filter(tib, as.logical(val)) %>% select(items) %>% rename(val = items)
        tibble(tag, val)
}

samples <- 
        tibble(files = list.files("./data/imvigor210/samples/"),
               path = paste0("./data/imvigor210/samples/", files),
               dat = map(path, parse_dat)) %>% 
        unnest(cols = c(dat)) %>% 
        pivot_wider(names_from = tag, values_from = val)

names <- 
        tibble(files = list.files("./data/imvigor210/samples/"),
               path = paste0("./data/imvigor210/samples/", files),
               dat = map(path, parse_dat_name)) %>% 
        unnest(cols = c(dat))

both <- full_join(samples, names, by = "files") %>% 
        select(-c(path.x, path.y, files, phenotype))


# Write -------------------------------------------------------------------

write_tsv(both, "./data/imvigor210/clin.tsv")
