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

(comp %>%
    mutate(GeneRatio = str_split(GeneRatio, '/', simplify = TRUE)[,2],
           GeneRatio = Count/as.numeric(GeneRatio),
           p.adjust = ifelse(p.adjust > .01, round(p.adjust, 2), '<0.01')) %>%
    ggplot(aes(x=Description, y = GeneRatio)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~ontology+Cluster, ncol = 1, scales = 'free_y', strip.position = 'right') +
    lims(y = c(-.2, .2)) +
    theme_bw() +
    geom_text(aes(y = -.05, label = Count)) +
    geom_text(aes(y = -.15, label = p.adjust)) +
    scale_y_continuous(limits = c(-.2, .2),
                       breaks = c(-.15, -.05, 0, .1, .2),
                       labels = c('FDR', '(n)', '0', '0.1\nFraction of\ngenes in module', 0.2),
                       name = '') +
    labs(x = '')) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'overrep_by_modules.png', sep = '/'),
         width = 20, height = 30, units = 'cm')
