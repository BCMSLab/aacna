# load libraries
library(tidyverse)
library(reshape2)
library(ggdendro)
library(cowplot)
library(ggraph)
library(tidygraph)
library(WGCNA)

# load data
load('data/wgcna.rds')

# global variables
figures_dir = 'manuscript/figures/'

# generate figures
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

(plot_grid(sft$fitIndices %>%
             mutate(y = -sign(slope)*SFT.R.sq) %>%
             ggplot(aes(x = Power, y = y)) +
             geom_point() +
             geom_hline(yintercept = 0.60172738, lty = 2, color = 'red') +
             labs(x = 'Soft Thresold (power)',
                  y = 'Scale-Free Topology') +
             theme_bw(),
           sft$fitIndices %>%
             ggplot(aes(x = Power, y = `mean.k.`)) +
             geom_point() +
             geom_vline(xintercept = 5, lty = 2, color = 'red') +
             labs(x = 'Soft Threshold (power)',
                  y = 'Mean Connectivity') +
             theme_bw(),
           labels = 'AUTO',
           label_size = 10,
           label_fontface = 'plain',
           nrow = 1,
           scale = .9)) %>%
  ggsave(plot = .,
         filename = paste(figures_dir, 'sft_power.png', sep = '/'),
         width = 16, height = 7, units = 'cm')

png(paste(figures_dir, 'gene_membership.png', sep = '/'),
    width = 10, height = 10, units = 'cm',
    res = 600, pointsize = 1)

net$diss %>%
  melt %>%
  mutate(value = ifelse(value == 0, NA, value)) %>%
  acast(Var1 ~ Var2) %>%
  TOMplot(dendro = hclust(as.dist(net$diss), method = 'average'),
          Colors = net$merged_colors)

dev.off()

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
