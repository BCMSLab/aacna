# load libraries
library(tidyverse)
library(xtable)
library(reshape2)

# load data
load('data/wgcna.rda')

# global variables
tables_dir <- 'manuscript/tables/'

# generate table
df <- interactions_evidence %>%
  setNames(c('to', 'from', 'evidence', 'value')) %>%
  bind_rows(interactions_evidence) %>%
  filter(evidence != 'combined_score') %>%
  mutate(evidence = str_split(evidence, '\\_', simplify = TRUE)[,1]) %>%
  select(-value) %>%
  unique()

ampk <- unique(split(ann$symbol, ann$category)$AMPK)

networks %>%
  mutate_all(as.character) %>%
  transform(from = ifelse(from %in% ampk, from, to),
            to = ifelse(from %in% ampk, to, from)) %>%
  filter(from %in% ampk) %>%
  left_join(gene$df, by = c('from'='symbol')) %>%
  left_join(gene$df, by = c('to'='symbol')) %>%
  mutate(color = ifelse(color.x == color.y, color.x, 'inbetween')) %>%
  left_join(df) %>%
  mutate(num = as.integer(as.factor(evidence))) %>%
  mutate(tnote = ifelse(num == 'NA', '', paste('\\tnotex{', num, '}', sep = ''))) %>%
  group_by(from, to) %>%
  mutate(notes = ifelse(is.na(tnote), '', paste(tnote, collapse = ' \\tnote{,} '))) %>%
  ungroup() %>%
  mutate(to = paste(to, notes, sep = '')) %>%
  ungroup() %>%
  with({
      tab <- group_by(., color, from) %>%
        summarise(to = paste(unique(to), collapse = ', ')) %>%
        ungroup() %>%
        mutate(color = ifelse(duplicated(color), '', color)) %>%
        setNames(c('Module/color', 'AMPK', 'Autophagy')) %>%
        xtable(caption = 'Summary of reported and novel AMPK-autophagy interactions.',
               label = 'tab:module_interactions',
               align = 'cllp{.6\\textwidth}')

      means <- data_frame(
        evidence = c('coexpression', 'database', 'experiments', 'textmining'),
        means = c('in the same or in other species (transferred by homology)',
                  'gathered from curated databases',
                  'gathered from other protein-protein interaction databases',
                  'extracted from the abstracts of scientific literature')
      )

      com <- select(., num, evidence) %>%
        unique() %>%
        na.omit() %>%
        arrange(num) %>%
        left_join(means) %>%
        with(
          paste('\\item[', num, '] \\label{', num, '} \\textit{', evidence, '} ', means, sep = '')) %>%
        paste(collapse = '\n')
      add.to.row <- list(pos = list(6),
                         command = paste("\\begin{tablenotes}",
                                         com,
                                         "\\end{tablenotes}\n",
                                         sep = '\n'))

      fl <- paste(tables_dir, 'module_interactions.tex', sep = '')

      tab %>%
        print(include.rownames = FALSE,
              floating.environment = 'threeparttable',
              add.to.row = add.to.row,
              booktabs = TRUE,
              caption.placement = 'top',
              sanitize.text.function = identity,
              file = fl)

      txt <- read_lines(fl)
      ind <- grep('tablenotes', txt)
      after <- txt[(ind[2]+1):(ind[2]+2)]
      notes <- txt[ind[1]:ind[2]]
      txt[ind[1]:(ind[1]+1)] <- after
      txt[(ind[1]+2):(ind[2]+2)] <- notes
      write_lines(txt, fl)
  })

