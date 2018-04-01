# load libraries
library(tidyverse)
library(reshape2)
library(ggdendro)
library(cowplot)
library(ggraph)
library(tidygraph)
library(WGCNA)

# load data
load('data/wgcna.rda')

# global variables
figures_dir = 'manuscript/figures/'

# generate figure

(module_preserv$quality$observed$ref.GSE34150[-1] %>%
    map(rownames_to_column, var = 'color') %>%
    bind_rows(.id = 'gse') %>%
    mutate(gse = str_split(gse, '\\.', simplify = TRUE)[, 2]) %>%
    ggplot(aes(x = moduleSize, y = medianRank.qual, color = color)) +
    geom_point() +
    facet_wrap(~gse, nrow = 1) +
    theme_bw() +
    scale_color_manual(values = c('blue', 'gold', 'grey', 'turquoise')) +
    theme(legend.position = 'top') +
    labs(x = 'Module Size', y = 'Preservation Median Rank', color = '') +
    guides(color = guide_legend(nrow = 1))) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'preservation_ranks.png', sep = '/'),
         width = 18, height = 9, units = 'cm',
         pointsize = 9)
