# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table
#
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
