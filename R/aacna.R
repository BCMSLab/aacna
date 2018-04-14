#' \code{aacna} package
#'
#' Analysis Compendium of the Autophagy and AMPK Co-expression
#'
#' @docType package
#' @name aacna
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
## fix by @jennybc
## source https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
if(getRversion() >= "2.15.1")  utils::globalVariables(c('.','symbol','GO','from', 'to', 'value'))
