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

set.seed(12345)

(mutate(interactions, Interaction = 'STRING') %>%
    right_join(networks) %>%
    mutate(Interaction = ifelse(is.na(Interaction), 'Novel', Interaction)) %>%
    as_tbl_graph() %>%
    left_join(gene$df, by = c('name' = 'symbol')) %>%
    left_join(dplyr::select(ann, symbol, category) %>%
                unique() %>%
                setNames(c('name', 'Gene'))) %>%
    filter(color != 'grey') %>%
    filter(!name %in% c('Atg5', 'Cln3', 'Zc3h12a')) %>%
    with(
      map(c('blue', 'turquoise'), function(x) {
        filter(., color == x) %>%
          ggraph(layout = 'kk') +
          geom_edge_link(aes(color = Interaction)) +
          geom_node_point(size = 5, aes(color = Gene)) +
          scale_color_manual(values = c('green', 'gray')) +
          scale_edge_color_manual(values = c('lightgray', 'red')) +
          geom_node_text(aes(label = name), color = 'royalblue') +
          theme_graph() +
          theme(legend.position = 'top',
                legend.direction = 'vertical')
      })
    ) %>%
    plot_grid(plotlist = .,
              labels = 'AUTO',
              label_size = 12,
              label_fontface = 'plain')) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'module_networks.png', sep = '/'),
         width = 30, height = 20, units = 'cm',
         pointsize = 7)
