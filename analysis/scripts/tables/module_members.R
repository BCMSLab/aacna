# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table

gene$df %>%
  left_join(ann) %>%
  group_by(color, category) %>%
  summarise(gene = paste(unique(symbol), collapse = ', ')) %>%
  spread(category, gene) %>%
  na.omit() %>%
  setNames(c('Module/color', 'AMPK', 'Autophagy')) %>%
  xtable(caption = 'AMPK and autophagy genes in different modules/colors.',
         label = 'tab:module_members',
         align = 'clp{.2\\textwidth}p{.6\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_members.tex', sep = '/'))
