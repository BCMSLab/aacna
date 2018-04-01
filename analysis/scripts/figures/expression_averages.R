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

melt(multi_data) %>%
  dplyr::select(-L2) %>%
  setNames(c('gsm', 'symbol', 'exp', 'gse')) %>%
  group_by(gse, symbol) %>%
  summarise(ave = mean(exp)) %>%
  left_join(unique(gene$df)) %>%
  mutate(ave = ifelse(ave < 0, 0, ave)) %>%
  dcast(color + symbol ~ gse, value.var = 'ave') %>%
  mutate_at(vars(3:5), function(x) log(x+1)) %>%
  with({
    map(colnames(.)[c(3,4,6)], function(x){
      c <- stats::cor(.['GSE34150'], .[x], use = 'complete')
      ggplot(.) +
        geom_point(aes_string(x = 'GSE34150', y = x, color = 'color')) +
        theme_bw() +
        theme(legend.position = 'non',
              plot.title = element_text(size = 8)) +
        lims(y = c(0,10)) +
        scale_color_manual(values = unique(.$color)) +
        ggtitle(label = paste("Pearson's Coeff: ", round(c, 2)))
    })
  }) %>%
  plot_grid(plotlist = .,
            nrow = 1) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'expression_averages.png', sep = '/'),
         width = 18, height = 6, units = 'cm',
         pointsize = 1)
