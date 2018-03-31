# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rds')

# global variables
tables_dir <- 'manuscript/tables/'


# generate tables

ann %>%
  group_by(category, term) %>%
  summarise(genes = paste(unique(symbol), collapse = ', ')) %>%
  ungroup() %>%
  mutate(category = ifelse(duplicated(category), '', category)) %>%
  setNames(c('Category', 'Term', 'Genes')) %>%
  xtable(caption = 'Gene members of the AMPK and autopahgy gene ontology terms.',
         align = 'clp{4cm}p{9cm}',
         label = 'tab:go_genes') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        file = paste(tables_dir, 'go_genes.tex', sep = '/'))

data_frame('Series ID' = c('GSE15018', 'GSE20696', 'GSE34150', 'GSE69313'),
           'Platform ID' = c('GPL6845', 'GPL1261', 'GPL6885', 'GPL6246'),
           Samples = c('54', '8', '24', '48'),
           Included = c('18', '8', '24', '12'),
           '(Contact, year)' = c('(Chin, 2009)', '(Mikkelsen, 2010)', '(Irmler, 2011)', '(Renbin, 2015)'),
           Reference = c("\\cite{ChinK2010Dataset:Cells}",
                         "\\cite{Mikkelsen2010ComparativeAdipogenesis}",
                         "\\cite{Horsch2015Dataset:Adipocytes}",
                         "\\cite{Zhang2015Dataset:Differentiation}")) %>%
  setNames(c('Series ID', 'Platform ID', 'Samples', 'Included', '(Contact, year)', 'Reference')) %>%
  xtable(caption = 'MDI-induced 3T3-L1 microarrays datasets.',
         align = 'ccccclc',
         label = 'tab:datasets') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'datasets.tex', sep = '/'))

gene$df %>%
  left_join(ann) %>%
  group_by(color, category) %>%
  summarise(gene = paste(unique(symbol), collapse = ', ')) %>%
  spread(category, gene) %>%
  na.omit() %>%
  setNames(c('Module/color', 'AMPK', 'Autophagy')) %>%
  xtable(caption = 'AMPK and autophagy genes in different modules/colors.',
         label = 'tab:module_members',
         align = 'clp{.2\\textwidth}p{.6\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_members.tex', sep = '/'))

ampk <- unique(split(ann$symbol, ann$category)$AMPK)

networks %>%
  mutate_all(as.character) %>%
  transform(from = ifelse(from %in% ampk, from, to),
            to = ifelse(from %in% ampk, to, from)) %>%
  filter(from %in% ampk) %>%
  left_join(gene$df, by = c('from'='symbol')) %>%
  left_join(gene$df, by = c('to'='symbol')) %>%
  mutate(color = ifelse(color.x == color.y, color.x, 'inbetween')) %>%
  group_by(color, from) %>%
  summarise(to = paste(unique(to), collapse = ', ')) %>%
  ungroup() %>%
  mutate(color = ifelse(duplicated(color), '', color)) %>%
  setNames(c('Module/color', 'AMPK', 'Autophagy')) %>%
  xtable(caption = 'Novel AMPK and autophagy interactions.',
         label = 'tab:module_interactions',
         align = 'cllp{.5\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_interactions.tex', sep = '/'))

list(GSE15018 = 1:18,
     GSE20696 = rep(c(-48, 0, 48, 168), 2),
     GSE34150 = rep(c(0, 48, 96, 144, 192, 240, 336, 432), 3),
     GSE69313 = rep(c(0, 6, 24, 72), 3)) %>%
  melt(value.name = 'Time Point (hours)') %>%
  dcast(`Time Point (hours)` ~ L1) %>%
  mutate_at(vars(-1), function(x) ifelse(x == 0, '', x)) %>%
  mutate(`Time Point (hours)` = as.character(`Time Point (hours)`)) %>%
  xtable(caption = 'Sample time points in the different datasets.',
         label = 'tab:time_points',
         align = 'cccccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        file = paste(tables_dir, 'time_points.tex', sep = '/'))

importance %>%
  group_by(module) %>%
  arrange(desc(hub)) %>%
  dplyr::slice(1:5) %>%
  ungroup() %>%
  mutate(module = ifelse(duplicated(module), '', module)) %>%
  setNames(c('Module/color', 'Gene', 'Degree', 'Betweenness', 'Closeness','Hub Score')) %>%
  xtable(caption = 'Top five hubs in the different module networks.',
         label = 'tab:network_hubs',
         align = 'cllcccc') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        file = paste(tables_dir, 'network_hubs.tex', sep = '/'))


df1 <- comp %>%
  transform(geneID = strsplit(geneID, '/')) %>%
  unnest(geneID)

df2 <- networks %>%
  mutate_all(as.character) %>%
  transform(from = ifelse(from %in% ampk, from, to),
            to = ifelse(from %in% ampk, to, from)) %>%
  filter(from %in% ampk) %>%
  left_join(gene$df, by = c('from'='symbol')) %>%
  left_join(gene$df, by = c('to'='symbol')) %>%
  mutate(color = ifelse(color.x == color.y, color.x, 'inbetween')) %>%
  dplyr::select(color, from, to)

inner_join(df1, df2, by =c('Cluster'='color', 'geneID'='to')) %>%
  group_by(Cluster, ontology, from, Description) %>%
  summarise(geneID = paste(geneID, collapse = ', ')) %>%
  ungroup() %>%
  group_by(Cluster) %>%
  mutate_at(vars(ontology, from), function(x) ifelse(duplicated(x), '', x)) %>%
  ungroup() %>%
  mutate(Cluster = ifelse(duplicated(Cluster), '', Cluster)) %>%
  setNames(c('Module/color', 'Ontology', 'AMPK', 'Term', 'Autophgy')) %>%
  xtable(caption = 'AMPK and autophagy interactions by gene ontology term.',
         label = 'tab:module_interactions_term',
         align = 'cllllp{.1\\textwidth}') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        table.placement = 'H',
        caption.placement = 'top',
        sanitize.text.function = identity,
        file = paste(tables_dir, 'module_interactions_term.tex', sep = '/'))
