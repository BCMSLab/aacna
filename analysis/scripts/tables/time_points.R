# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table

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
