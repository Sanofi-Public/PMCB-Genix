# Documentation and definitions for data and constants

#### Data ####

#' Example for demo
#'
#' A small example database subset from Tabula Sapiens Consortium, 2022.
#'
#' @format   An object with the following items:
#'   \itemize{
#'     \item  \code{gex}:               numeric matrix of single-cell 
#'                                      expression values where rows are 
#'                                      genes and columns are cells.
#'     \item  \code{mt.data}:           data.frame of meta data including
#'                                      organ tissue (\code{organ_tissue}),
#'                                      expert annotation (\code{free_annotation}), 
#'                                      and cell-type compartment
#'                                      abbreviations (\code{celltype}): "B" 
#'                                      for B cells, "MPh" for Myloides, 
#'                                      "TNK" for T and Natural killer cells, 
#'                                      "Endo" for Endothelial, and 
#'                                      "Epit" for Epithelial.
#'     \item  \code{gene.sets}:         list of gene-markers. Each element of 
#'                                      list is named by the cell-type 
#'                                      and contains a vector of associated 
#'                                      gene-markers.
#'     \item  \code{adjacent.mtx}:      numeric matrix representing the 
#'                                      nearest-neighbor graph. Matrix is
#'                                      square with number of rows and coulmns 
#'                                      equal to the number of cells.
#'     \item  \code{gene.tnk.sets}:     similar to \code{gene.sets}, but 
#'                                      specific to TNK subset
#'     \item  \code{adjacent.tnk.mtx}:  similar to \code{adjacent.mtx}, but 
#'                                      specific to TNK subset.
#' }
#' 
#' @return An R object. 
#' 
#' @references
#' \enumerate{
#'   \item  Tabula Sapiens Consortium, et al. The Tabula Sapiens: 
#'            A multiple-organ, single-cell transcriptomic atlas of humans.
#'            Science 376.6594 (2022): eabl4896.
#' }
#' 
#' @usage  data("genixDemoData") 
"genixDemoData"