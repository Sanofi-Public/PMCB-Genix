# Documentation and definitions for data and constants

#### Data ####

#' Example for demo
#'
#' A small example database containing the gene expression matrix of 200 cells.
#'
#' @format   An object with the following items:
#'   \itemize{
#'     \item  \code{gex}:               numeric matrix of single-cell 
#'                                      expression values where rows are 
#'                                      genes and columns are cells.
#'     \item  \code{metadata}:          data.frame of meta data including 
#'                                      typical \code{Seurat} metadata
#'                                      columns such as \code{orig.ident},
#'                                      \code{nCount_RNA}, and \code{nFeature_RNA}.
#'                                      Additionally, it includes a \code{group}
#'                                      column that identifies the group to
#'                                      which each cell belongs.
#' }
#' 
#' @return An R object. 
#' 
#' @usage  data("genixDemoData") 
"genixDemoData"