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

plot_grid(cor(net$mes, as.numeric(stage), use = 'p') %>%
            as.data.frame() %>%
            rownames_to_column('color') %>%
            setNames(c('color', 'cor')) %>%
            mutate(color = c('blue', 'turquoise', 'gray')) %>%
            ggplot(aes(x = color, y = cor)) +
            geom_col() +
            theme_bw() +
            lims(y = c(-1,1)) +
            labs(x = '', y = "Pearson's Correlation") +
            geom_abline(intercept = 0, slope = 0, lty = 2),
          mr %>%
            rownames_to_column('module') %>%
            mutate(module = c('blue', 'gray', 'turquoise')) %>%
            gather(dir, prop, starts_with('Prop')) %>%
            mutate(prop = ifelse(dir == 'PropDown', -prop, prop)) %>%
            ggplot(aes(x = module, y = prop)) +
            geom_col() +
            lims(y = c(-1,1)) +
            labs(x = '', y = 'Fraction of DE genes') +
            theme_bw() +
            theme(legend.position = 'top') +
            geom_abline(intercept = 0, slope = 0, lty = 2),
          labels = 'AUTO',
          label_fontface = 'plain',
          label_size = 10,
          scale = .9) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'module_cor_overrep.png', sep = '/'),
         width = 16, height = 7, units = 'cm')
