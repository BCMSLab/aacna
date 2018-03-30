context("functions")

test_that('test annotation_get', {
  library(org.Mm.eg.db)
  library(GO.db)

  go_ids <- c('GO:0006914', 'GO:0004679')
  go_names <- c('autophagy', 'ampk')
  go_db <- GOBPCHILDREN
  go_term <- GOTERM
  org_db <- org.Mm.eg.db
  columns <- 'SYMBOL'

  df <- annotation_get(go_ids,
                       go_names,
                       go_db,
                       go_term,
                       org_db,
                       columns)

  expect_s3_class(df, 'data.frame')
  expect_identical(unique(df$category), go_names)
  expect_true(all(!is.na(df$category)))
})

test_that('test expression_subset', {
  library(Biobase)
  data("sample.ExpressionSet")
  pd <- pData(sample.ExpressionSet)

  set.seed(123)
  byrow <- sample(c(TRUE, FALSE), nrow(sample.ExpressionSet), replace = TRUE)
  bycol <- pd$sex == 'Female'

  mat <- expression_subset(sample.ExpressionSet, byrow, bycol)

  expect_equal(dim(mat), c(sum(byrow), sum(bycol)))

  # test collapsing rows
  set.seed(123)
  rowGroup <- sample(LETTERS, nrow(sample.ExpressionSet), replace = TRUE)

  mat <- expression_subset(sample.ExpressionSet,
                           collapse_rows = TRUE, rowGroup = rowGroup,
                           rowID = featureNames(sample.ExpressionSet))

  expect_equal(dim(mat), c(ncol(sample.ExpressionSet)[[1]], length(unique(rowGroup))))
})

test_that('test expression_reshape', {
  library(Biobase)
  data("sample.ExpressionSet")

  mat <- exprs(sample.ExpressionSet)

  dat <- expression_reshape(mat)

  expect_true(all(dim(dat$data) == dim(mat)[2:1]))
})

test_that('test metadata_get', {
  library(Biobase)
  data("sample.ExpressionSet")

  df <- metadata_get(sample.ExpressionSet, c('sex', 'type'), c('sex', 'type'))

  expect_equal(dim(df), c(ncol(sample.ExpressionSet)[[1]], 3))
})


test_that('test cna_run', {
  library(Biobase)
  data("sample.ExpressionSet")
  dat <- list(data = t(exprs(sample.ExpressionSet)))

  cna <- cna_run(dat$data, power = 5)

  expect_equal(length(cna), 9)
  expect_equal(length(cna$colors), nrow(sample.ExpressionSet)[[1]])
})

test_that('test module_compare', {
  library(org.Mm.eg.db)
  set.seed(1234)
  genes <- sample(keys(org.Mm.eg.db, 'SYMBOL'), 10)
  colors <- sample(c('red', 'blue'), 10, replace = TRUE)
  ind <- split(genes, colors)

  comp <- module_compare(ind,
                         fun = 'enrichGO',
                         OrgDb = org.Mm.eg.db,
                         keyType = 'SYMBOL',
                         ont = 'CC',
                         pAdjustMethod = 'fdr')

  expect_s3_class(comp, 'data.frame')
  expect_identical(levels(comp$Cluster), c('blue', 'red'))
  })
