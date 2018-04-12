# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table

ampk <- unique(split(ann$symbol, ann$category)$AMPK)

df1 <- comp %>%
  transform(geneID = strsplit(geneID, '/')) %>%
  unnest(geneID)

df2 <- networks %>%
  mutate_all(as.character) %>%
  transform(from = ifelse(from %in% ampk, from, to),
            to = ifelse(from %in% ampk, to, from)) %>%
  filter(from %in% ampk) %>%
  left_join(gene$df, by = c('from'='symbol')) %>%
  left_join(gene$df, by = c('to'='symbol')) %>%
  mutate(color = ifelse(color.x == color.y, color.x, 'inbetween')) %>%
  dplyr::select(color, from, to)

inner_join(df1, df2, by =c('Cluster'='color', 'geneID'='to')) %>%
  group_by(Cluster, ontology, from, Description) %>%
  summarise(geneID = paste(geneID, collapse = ', ')) %>%
  ungroup() %>%
  group_by(Cluster) %>%
  mutate_at(vars(ontology, from), function(x) ifelse(duplicated(x), '', x)) %>%
  ungroup() %>%
  mutate(Cluster = ifelse(duplicated(Cluster), '', Cluster)) %>%
  setNames(c('Module/color', 'Ontology', 'AMPK', 'Term', 'Autophgy')) %>%
  xtable(caption = 'AMPK and autophagy interactions by gene ontology term.',
         label = 'tab:module_interactions_term',
         align = 'clllp{.4\\textwidth}p{.15\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_interactions_term.tex', sep = '/'))
