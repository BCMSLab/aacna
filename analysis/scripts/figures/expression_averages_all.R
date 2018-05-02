# load libraries
library(tidyverse)
library(reshape2)
library(cowplot)
library(Biobase)
library(aacna)
library(GEOquery)

# load data
load('data/wgcna.rda')

# global variables
figures_dir = 'manuscript/figures/'

# generate figure
esets <- map2(gse, samples, function(x, y) {
  if(x == 'GSE69313') {
    e <- getGEO('GSE69313', destdir = data)[[1]]

    mo <- mogene10sttranscriptcluster.db::mogene10sttranscriptcluster.db
    GPL6246 <- AnnotationDbi::select(mo, featureNames(e),
                                     'SYMBOL', 'PROBEID',
                                     multiVals=first) %>%
      dplyr::filter(!duplicated(PROBEID))
    all(featureNames(e) == GPL6246$PROBEID)
    rownames(GPL6246) <- GPL6246$PROBEID
    fData(e) <- GPL6246

    eset <- e
  } else {
    eset <- GEOquery::getGEO(x, destdir = data, AnnotGPL = FALSE)[[1]]
  }


  col <- 1:ncol(eset) %in% y

  col_name <- grep('symbol', names(fData(eset)), ignore.case = TRUE, value = TRUE)
  fd <- metadata_get(eset, col_name, 'symbol', type = 'features')

  e <- expression_subset(eset,
                         collapse_rows = TRUE,
                         rowGroup = fd$symbol, rowID = fd$probe_id)
  return(e)
})

names(esets) <- gse

esets %>%
  map(melt) %>%
  bind_rows(.id = 'gse') %>%
  setNames(c('gse', 'symbol', 'gsm', 'exp')) %>%
  group_by(gse, symbol) %>%
  summarise(ave = mean(exp)) %>%
  mutate(ave = ifelse(ave < 0, 0, ave)) %>%
  dcast(symbol ~ gse, value.var = 'ave') %>%
  mutate_at(vars(2:4), function(x) log(x+1)) %>%
  with({
    map(colnames(.)[c(2, 3, 5)], function(x){
      c <- stats::cor(.['GSE34150'], .[x], use = 'complete')
      ggplot(.) +
        geom_point(aes_string(x = 'GSE34150', y = x), alpha = .2, na.rm = TRUE) +
        theme_bw() +
        theme(legend.position = 'non',
              plot.title = element_text(size = 8)) +
        lims(y = c(0,10)) +
        ggtitle(label = paste("Pearson's Coeff: ", round(c, 2)))
    })
  }) %>%
  plot_grid(plotlist = .,
            nrow = 1) %>%
ggsave(plot = .,
       filename = paste(figures_dir, 'expression_averages_all.png', sep = '/'),
       width = 18, height = 6, units = 'cm',
       pointsize = 1)

