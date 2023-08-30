# Classes, Generics and Methods

#### Generics ####

setOldClass("igraph")
setClassUnion("ListIgraph", members=c("list", "igraph"))


#### Genix classes ####

#' S4 class defining genixObj class
#'
#' \code{genixObj} defines constructNets outputs.
#' 
#' @slot  eic_graph   an estimated inverse covariance graph. The edges indicate
#'                    the partial correlation values. The nodes represent the genes.
#' 
#' @name         genixObj-class
#' @rdname       genixObj-class
#' @aliases      genixObj
#' @exportClass  genixObj
#' @return       genixObj
setClass("genixObj", 
         slots=c(icm_graph = "ListIgraph"))


#' S4 class defining compiledObj class
#'
#' \code{compiledObj} defines compileNets outputs.
#' 
#' @slot  graph             a compiled estimated inverse covariance graph. The edges indicate
#'                          the partial correlation values. The nodes represent the genes.
#' @slot  graph_df          converted graph to a long data frame:
#'  \itemize{ 
#'            \item    \code{from_gene}: vertex A. 
#'            \item    \code{from_connectivity}: vertex A degree. 
#'            \item    \code{from_module}: the module vertex A belongs to. 
#'            \item    \code{from_hub}: is the vertex A, a hub gene. 
#'            \item    \code{par_corr_coeff}: The weight of given edge equal to 
#'                                            the partial correlation coefficient. 
#'            \item    \code{to_gene}: vertex B. 
#'            \item    \code{to_connectivity}: vertex B degree. 
#'            \item    \code{to_module}: the module vertex B belongs to. 
#'            \item    \code{to_hub}: is the vertex B, a hub gene. 
#'            }  
#' @slot  hubs               A named vector containing the names of hub genes
#'                           and their connectivity, sorted in ascending order.
#' @slot  modules            A named list containing the module gene names and  
#'                           their connectivity, sorted from lowest degree 
#'                           ascending. The modules themselves are sorted in descending 
#'                           order by size, starting from the largest module (M0).
#' 
#' @name         compiledObj-class
#' @rdname       compiledObj-class
#' @aliases      compiledObj
#' @exportClass  compiledObj
#' @return       compiledObj
setClass("compiledObj", 
         slots=c(graph = "igraph",
                 graph_df = "data.frame",
                 hubs = "numeric",
                 modules = "list"))


#' S4 class defining comparedObj class
#'
#' \code{comparedObj} defines compareNets outputs.
#' 
#' @slot  ji             The Jaccard index between grph.1 and grph.2,
#'                       considering both shared genes and edges.
#' @slot permutations    A numeric vector of calculated Jaccard indices,
#'                       comparing pairs of randomly reconstructed networks that 
#'                       are sampled from a union of observed grph.1 and grph.2.
#' @slot pvalue          The \code{permutations} distribution is used to calculate a 
#'                       two-sided t-test p-value, which reflects the proportion 
#'                       of permuted Jaccard indices that differ from the observed \code{ji}.
#' @slot  deltas_df      A data frame containing the features calculated per removal 
#'                       of a given gene by comparing two networks grph.1 and grph.2.
#'  \itemize{ 
#'            \item    \code{gene}: gene name. 
#'            \item    \code{jaccard}: The Jaccard similarity index between grph.1 and grph.2, 
#'            after gene removal.
#'            \item    \code{degree.1} and \code{degree.2}: the normalized degree of a
#'            given gene in grph.1 and grph.2. Zero, if the gene does not exist.
#'            \item    \code{closeness.1} and \code{closeness.2}: the normalized closeness
#'            of a given gene in grph.1 and grph.2. Zero, if the gene does not
#'            exist.
#'            \item    \code{betweenness.1} and \code{betweenness.2}: the normalized betweenness
#'            of a given gene in grph.1 and grph.2. Zero, if the gene does not exist.
#'            \item    \code{delta_jaccard}: the relative difference in similarity of 
#'            grph.1 and grph.2 after removal of given gene (\code{jaccard}) compared 
#'            to before removal (\code{ji}).
#'            \item    \code{delta_degree}: the difference between degree.1 and degree.2.
#'            \item    \code{delta_closeness}: the difference between closeness.1 and closeness.2
#'            \item    \code{delta_betweenness}: the difference between betweenness.1 and betweenness.2.
#'            \item    \code{type}: wether a gene is hub in grph.1, grph.2, both, or non.
#'            }  
#' @slot  inA             genes in grph.1 and not in grph.2
#' @slot  inB             genes in grph.2 and not in grph.1
#' @slot  inAB            genes in both grph.1 and grph.2
#' @slot  inAB_mtx        table of gene counts in grph.1 and grph.2
#'  
#' @name         comparedObj-class
#' @rdname       comparedObj-class
#' @aliases      comparedObj
#' @exportClass  comparedObj
#' @return       comparedObj
setClass("comparedObj", 
         slots=c(ji = "numeric",
                 pvalue = "numeric",
                 permutatins = "numeric",
                 deltas = "data.frame",
                 inA = "character",
                 inB = "character",
                 inAB = "character",
                 inAB_mtx = "table"))
