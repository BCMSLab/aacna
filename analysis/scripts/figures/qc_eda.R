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
(plot_grid(
  melt(mat) %>%
    ggplot(aes(x=Var2, y=value)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = .5, size = 8)) +
    labs(x = '', y = 'Log expression'),
  mat %>%
    t() %>%
    dist() %>%
    hclust() %>%
    as.dendrogram() %>%
    dendro_data() %>%
    with({
      df = .$segment %>%
        mutate(yend = ifelse(yend == 0, 80, yend))

      pd <- data_frame(stage = stage, label = colnames(mat))

      labels = .$labels %>%
        left_join(pd)

      df %>%
        ggplot() +
        geom_segment(aes(x=x,y=y,xend=xend,yend=yend)) +
        geom_text(data = labels,
                  aes(x=x,y=75, label=as.character(label), color=stage),
                  angle = 90, hjust = 1, size = 3) +
        lims(y=c(58, 180)) +
        geom_abline(slope = 0, intercept = 120, lty=2) +
        theme_dendro() +
        theme(legend.position = 'none',
              plot.margin = margin(-2,0,0,0, 'cm'))
    }),
  cmdscale(dist(t(mat))) %>%
    as.data.frame() %>%
    ggplot(aes(x=V1, y=V2, color = stage)) +
    geom_point() +
    theme_bw() +
    labs(x = 'D1', y = 'D2') +
    theme(legend.position = 'none'),
  labels = 'AUTO',
  label_size = 10,
  label_fontface = 'plain',
  nrow = 1,
  scale = .9
)) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'qc_eda.png', sep = '/'),
         width = 30, height = 10, units = 'cm')
