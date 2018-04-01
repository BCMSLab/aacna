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

(plot_grid(similarity %>%
             filter(abs(value) < 1) %>%
             ggplot(aes(x = abs(value), color = measure)) +
             stat_ecdf() +
             theme_bw() +
             theme(legend.position = 'top',
                   legend.direction = 'vertical') +
             labs(y = 'CDF',
                  x = 'Correlation',
                  color = 'Measure'),
           importance %>%
             ggplot(aes(x=degree, y=betweenness, size = hub)) +
             geom_point() +
             labs(x = 'Degree Centrality',
                  y = 'Betweenness Centrality',
                  size = 'Hub Score') +
             scale_size_continuous(breaks = c(0, .5, 1)) +
             theme_bw() +
             theme(legend.position = 'top', legend.direction = 'vertical'),
           labels = 'AUTO',
           label_size = 10,
           label_fontface = 'plain',
           nrow = 1,
           scale = .9)) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'similarity_centrality.png', sep = '/'),
         width = 18, height = 12, units = 'cm')
