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

(module_preserv$quality$Z$ref.GSE34150[-1] %>%
    map(rownames_to_column, var = 'color') %>%
    bind_rows(.id = 'gse') %>%
    mutate(gse = str_split(gse, '\\.', simplify = TRUE)[, 2]) %>%
    ggplot(aes(x = moduleSize, y = Zsummary.qual, color = color)) +
    geom_point() +
    geom_abline(intercept = c(2,5), slope = 0, lty = 2) +
    theme_bw() +
    scale_color_manual(values = c('blue', 'gold', 'grey', 'turquoise')) +
    theme(legend.position = 'top') +
    facet_wrap(~gse, nrow = 1) +
    labs(x = 'Module Size', y = 'Preservation Z Summary', color = '') +
    guides(color = guide_legend(nrow = 1))) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'preservation_zsummary.png', sep = '/'),
         width = 18, height = 9, units = 'cm')
