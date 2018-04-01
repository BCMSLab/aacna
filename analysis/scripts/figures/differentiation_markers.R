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

plot_grid(
  mat %>%
    melt() %>%
    filter(Var1 %in% c('Lpl', 'Cebpa', 'Pparg')) %>%
    left_join(data_frame(stage = stage,
                         Var2 = colnames(mat))) %>%
    group_by(stage, Var1) %>%
    summarise(ave = mean(value)) %>%
    ggplot(aes(x = Var1, y = ave, group = stage, fill = stage)) +
    geom_col(position = 'dodge') +
    theme_bw() +
    theme(legend.position = 'top') +
    labs(x = '',
         y = 'Average log expression',
         fill = '') +
    guides(fill = guide_legend(nrow = 3)),
  mat %>%
    melt() %>%
    filter(Var1 %in% c('Fasn', 'Acly', 'Acaca', 'Elovl6', 'Scd1', 'Scd2',
                       'Scd3', 'Scd4', 'Dgat1', 'Dgat2')) %>%
    left_join(data_frame(day = as.factor(rep(c(seq(0, 10, 2), 14, 18), 3)),
                         Var2 = colnames(mat))) %>%
    group_by(day, Var1) %>%
    summarise(ave = mean(value)) %>%
    ggplot(aes(x = day, y = ave, group = Var1, color = Var1)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = 'top') +
    labs(x = 'Time point (Day)',
         y = 'Average log expression',
         color = '') +
    guides(color = guide_legend(nrow = 3)),
  labels = 'AUTO',
  label_fontface = 'plain',
  label_size = 10,
  scale = .9
) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'differentiation_markers.png', sep = '/'),
         width = 18,
         height = 12,
         units = 'cm')
