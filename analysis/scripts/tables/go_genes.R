# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table

ann %>%
  group_by(category, term) %>%
  summarise(genes = paste(unique(symbol), collapse = ', ')) %>%
  ungroup() %>%
  mutate(category = ifelse(duplicated(category), '', category)) %>%
  setNames(c('Category', 'Term', 'Genes')) %>%
  xtable(caption = 'Gene members of the AMPK and autopahgy gene ontology terms.',
         align = 'clp{4cm}p{9cm}',
         label = 'tab:go_genes') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        file = paste(tables_dir, 'go_genes.tex', sep = '/'))
