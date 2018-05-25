# load libraries
library(tidyverse)
library(readxl)
library(xtable)

# global variable
tables_dir <- 'manuscript/tables/'

# load data
primers <- read_excel('data/rtqpcr.xlsx', sheet = 3)

# generate tables
primers %>%
  xtable(caption = 'RT-qPCR primer sequences.',
         align = 'clll',
         label = 'tab:rtqpcr_primers') %>%
  print(include.rownames = FALSE,
        booktabs = TRUE,
        caption.placement = 'top',
        table.placement = 'H',
        file = paste(tables_dir, 'rtqpcr_primers.tex', sep = '/'))
