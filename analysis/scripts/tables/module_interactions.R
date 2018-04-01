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

networks %>%
  mutate_all(as.character) %>%
  transform(from = ifelse(from %in% ampk, from, to),
            to = ifelse(from %in% ampk, to, from)) %>%
  filter(from %in% ampk) %>%
  left_join(gene$df, by = c('from'='symbol')) %>%
  left_join(gene$df, by = c('to'='symbol')) %>%
  mutate(color = ifelse(color.x == color.y, color.x, 'inbetween')) %>%
  group_by(color, from) %>%
  summarise(to = paste(unique(to), collapse = ', ')) %>%
  ungroup() %>%
  mutate(color = ifelse(duplicated(color), '', color)) %>%
  setNames(c('Module/color', 'AMPK', 'Autophagy')) %>%
  xtable(caption = 'Novel AMPK and autophagy interactions.',
         label = 'tab:module_interactions',
         align = 'cllp{.5\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_interactions.tex', sep = '/'))
