# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table

importance %>%
  group_by(module) %>%
  arrange(desc(hub)) %>%
  dplyr::slice(1:5) %>%
  ungroup() %>%
  mutate(module = ifelse(duplicated(module), '', module)) %>%
  setNames(c('Module/color', 'Gene', 'Degree', 'Betweenness', 'Closeness','Hub Score')) %>%
  xtable(caption = 'Top five hubs in the different module networks.',
         label = 'tab:network_hubs',
         align = 'cllcccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        file = paste(tables_dir, 'network_hubs.tex', sep = '/'))
