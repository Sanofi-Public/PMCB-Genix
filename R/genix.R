# Genix package documentation and import directives

#' The genix package
#' 
#' \code{genix} provides a foundation for sinature discovery through
#' comparative analysis of gene association networks using single cell RNA
#' sequencing data. For additional details regarding the use of the \code{genix}
#' package see the vignettes:\cr
#' \code{browseVignettes("genix")}
#' 
#' @section  Construction, Process, and Comparison of Gene Association Networks:
#' \itemize{
#'   \item  \link{constructNets}:     Construct gene association network from
#'                                   scRNA-seq data.
#'   \item  \link{compileNets}:       Compile the constructed network and
#'                                   identify hub gene and gene modules.
#'   \item  \link{compareNets}:       Compare two ntworks and report topological 
#'                                   variant genes.
#' }
#' 
#' @name     genix
#' @docType  package
#' @references
#' \enumerate{
#'   \item  Nima Nouri, et al.: Signature Genes Identification through
#'                              Comparative Analysis of Gene Association
#'                              Networks Using Single-Cell RNA Sequencing Data
#' }
#'
#' @import      ggplot2
#' @import      methods
#' @importFrom  glasso          glasso
#' @importFrom  dplyr           n %>% filter select mutate summarize rename
#' @importFrom  igraph          intersection union ecount vcount degree closeness 
#'                              betweenness transitivity graph_from_adjacency_matrix 
#'                              as_adjacency_matrix delete.vertices V E 
#'                              subgraph.edges
#' @importFrom  gtools          mixedorder
#' @importFrom  scales          log10_trans trans_breaks trans_format
#' @importFrom  Matrix          rowSums colSums Matrix
#' @importFrom  rlang           sym syms 
#' @importFrom  stats           t.test hclust as.dist cov IQR quantile
#' @importFrom  graphics        hist
#' @importFrom  dynamicTreeCut  cutreeDynamicTree
NULL

# Package loading actions
.onAttach <- function(libname, pkgname) {
    msg <- paste("genix package.", 
                 "License and Copyrights in file LICENSE and COPYRIGHTS.",
                 sep="\n")
    packageStartupMessage(msg)
}
