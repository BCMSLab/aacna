# load libraries
library(aacna)
library(WGCNA)
library(GO.db)
library(org.Mm.eg.db)
library(tidyverse)
library(GEOquery)
library(testthat)
library(limma)
library(igraph)

# global variables
data <- 'data/'
output <- 'wgcna.rds'

allowWGCNAThreads(2)

# get data
## annotation
go_ids <- c('GO:0006914', 'GO:0004679')
go_names <- c('autophagy', 'AMPK')
go_db <- GOBPCHILDREN
go_term <- GOTERM
org_db <- org.Mm.eg.db
columns <- 'SYMBOL'

ann <- annotation_get(go_ids,
                     go_names,
                     go_db,
                     go_term,
                     org_db,
                     columns)
# check
test_that('correct annoation', {
  expect_identical(unique(ann$category), go_names)

  ampk <- dplyr::filter(ann, category == 'AMPK') %>%
    pull(symbol) %>%
    unique()

  expect_true(length(ampk) == 14)

  autophagy <- dplyr::filter(ann, category == 'autophagy') %>%
    pull(symbol) %>%
    unique()

  expect_true(length(autophagy) == 167)
  })

## expression data
gse <- c('GSE34150', 'GSE15018', 'GSE20696', 'GSE69313')

fls <- paste(data, gse, '_series_matrix.txt.gz', sep = '')

if(all(!file.exists(fls))) {
  expression_get(gse,
                 destdir = data,
                 AnnotGPL = FALSE,
                 getGPL = FALSE)
}

# check
test_that('gse matrices exist', {
  expect_true(
    all(
      file.exists(fls)
    )
  )
})

# reshape data
## subset to relevant samples
samples <- list(1:24,
                1:18,
                1:8,
                c(1:3, 13:15, 25:27, 37:39))

genes <- unique(ann$symbol)

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
  row <- fd$symbol %in% genes

  e <- expression_subset(eset, row, col,
                         collapse_rows = TRUE,
                         rowGroup = fd$symbol, rowID = fd$probe_id)
  return(e)
})
names(esets) <- gse

# check esets
test_that('esets are subseted correctly', {
  expect_true(length(esets) == 4)
  expect_true(all(lengths(samples) == sapply(esets, dim)[2,]))
  expect_true(all(sapply(esets, rownames)[[4]] %in% ann$symbol))
})

# reshape data
multi_data <- map(esets, expression_reshape)

# check
test_that('multi_data in appropriate shape', {
  expect_true(
    all(
      map(multi_data, function(x) dim(x$data))[[1]] == map(esets, function(x) dim(x)[2:1])[[1]]
    )
  )
  expect_true(
    all(unlist(map(multi_data, function(x) sum(is.na(x$data)))) == 0)
  )
})

# check samples and genes
goodSamplesGenes(multi_data$GSE34150$data, verbose = FALSE)$allOK

# calculate threshold
sft <- pickSoftThreshold(multi_data$GSE34150$data)

# run analysis
net <- cna_run(multi_data$GSE34150$data, 5)

# check
test_that('cna ran successfully', {
  expect_equal(length(unique(net$merged_colors)), 3)
  expect_equal(length(net$merged_colors), ncol(multi_data$GSE34150$data))
})

# unpack members
gene <- list(
  df = data.frame(
    symbol = colnames(multi_data$GSE34150$data),
    color = net$merged_colors,
    stringsAsFactors = FALSE
  ),
  lst = split(colnames(multi_data$GSE34150$data), net$merged_colors)
)

# module correlation to stage

# module over-representation
e <- getGEO('GSE34150', destdir = data)[[1]]
mat <- collapseRows(exprs(e),
                    rowGroup = fData(e)$Symbol,
                    rowID = featureNames(e))[[1]]
mat[mat < 0] <- 0
mat <- log(mat + 1)

stage <- c('undifferentiated', 'differentiating', 'maturating')
stage <- rep(rep(stage, times = c(1,2,5)), 3)
stage <- factor(stage, levels = c('undifferentiated', 'differentiating', 'maturating'))
mod <- model.matrix(~stage)

mr <- mroast(mat,
             index = gene$lst,
             design = mod)

# check
test_that('overrepresentation done correctly', {
  expect_true(all(rownames(mr) %in% names(gene$lst)))
})

# check
# module to network
# edge list
networks <- exportNetworkToVisANT(net$adj, threshold = .1)

# string interactions
interactions <- interactions_get(gene$df,
                                 species = 10090,
                                 input_directory = tempdir())


test_that('string interactions obtained correctly', {
  expect_s3_class(interactions, 'data.frame')
  expect_identical(names(interactions), c('from', 'to'))
})

# make graph objects
gs <- module_network(gene$lst[-2], networks[, 1:2])

# check
test_that('network graphs constructed properly', {
  expect_identical(class(gs[[1]]), 'igraph')
  expect_true(all(igraph::V(gs$blue)$name %in% gene$lst$blue))
  expect_true(all(igraph::V(gs$turquoise)$name %in% gene$lst$turquoise))
  })

# member importance
importance <- member_importance(gs)

# pairwise similarity
similarity <- member_similarity(multi_data$GSE34150$data)

# module comparisons
comp <- list()

comp$MF <- module_compare(gene$lst[-2],
                          fun = 'enrichGO',
                          OrgDb = org.Mm.eg.db,
                          keyType = 'SYMBOL',
                          ont = 'MF',
                          pAdjustMethod = 'fdr')

comp$CC <- module_compare(gene$lst[-2],
                          level = 3,
                          fun = 'enrichGO',
                          OrgDb = org.Mm.eg.db,
                          keyType = 'SYMBOL',
                          ont = 'CC',
                          pAdjustMethod = 'fdr')

comp <- bind_rows(comp, .id = 'ontology')

# check
test_that('comarison done correctly', {
  expect_identical(unique(comp$ontology), c('MF', 'CC'))
  expect_identical(levels(comp$Cluster), names(gene$lst[-2]))
})

# module preservation
multi_color <- net$merged_colors
names(multi_color) <- colnames(multi_data$GSE34150$data)
multi_color <- list(GSE34150 = multi_color)
allowWGCNAThreads(3)
module_preserv <- modulePreservation(multiData = multi_data,
                                     multiColor = multi_color,
                                     nPermutations = 10)

# clean

# save
save.image('data/wgcna.rds')
