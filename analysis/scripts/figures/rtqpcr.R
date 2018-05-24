# load libraries
library(tidyverse)
library(reshape2)
library(readxl)
library(pcr)
library(cowplot)

# load data
qpcr <- map(list(b = 1, g = 2),
            function(x) read_excel('data/rtqpcr.xlsx', sheet = x))
norm_rel <- map(qpcr, function(x) {
  ct <- select(x, 3:ncol(x))
  group <- x$stage
  group <- factor(group, levels = unique(group))
  pcr_analyze(ct,
              group = group,
              reference_gene = '18s',
              reference_group = 'Confluent')
})

c(norm_rel %>%
  map(function(x) {
    x %>%
      ggplot(aes(x = group, y = relative_expression, group = gene, color = gene)) +
      geom_line() +
      geom_point() +
      theme_bw() +
      theme(legend.position = 'top',
            legend.direction = 'horizontal',
            axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(x = '',
           y = 'Relative mRNA level\n',
           color = '') +
      lims(y = c(0,3))
    }),
  norm_rel %>%
    map(function(x) {
      x %>%
        acast(group ~ gene, value.var = 'relative_expression') %>%
        cor() %>%
        melt() %>%
        filter(grepl('Prka*', Var1) & value != 1) %>%
        ggplot(aes(x = Var2, y = value)) +
        geom_col() +
        geom_hline(yintercept = 0, lty = 'dashed') +
        theme_bw() +
        labs(x = '',
             y = "Pearson's Coeffecient") +
        lims(y = c(-1, 1))
    })) %>%
  plot_grid(plotlist = .,
            labels = 'AUTO',
            label_fontface = 'plain',
            label_size = 11,
            scale = .9) %>%
  ggsave(filename = 'manuscript/figures/rtqpcr.png',
         width = 20, height = 20, units = 'cm')
