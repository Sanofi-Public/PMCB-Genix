# Main function and helper functions for genix

#' @include Classes.R
NULL


#### Main function ####

#' @title Construct gene association network
#' @description   A glasso-based algorithm to infer undirected dependency graphs 
#'                from scRNA-seq data
#'
#' @param    gex               A numeric matrix of single-cell expression values
#'                             where rows are genes and columns are cells. 
#' @param    glasso.rho        Regularization parameter for lasso. This scalar value 
#'                             controls the degree of sparsity in the resulting
#'                             graph. Higher values result in greater sparsity. 
#'                             See \link[glasso]{glasso}. Default is 1.
#' @param    metadata          Provide to split the input \code{gex} matrix by 
#'                             \code{field}. Row names should represent cell barcodes 
#'                             similar to the column names of the \code{gex} matrix. 
#' @param    field             The column name from \code{metadata} by which to
#'                             split the input \code{gex} matrix. This parameter 
#'                             results in one graph per group. Default is no spliting.
#' @param    scale.gex         if TRUE, the count expression matrix is normalized 
#'                             to the mean library size, where each cell is scaled 
#'                             to sum up to the mean total counts. Default is false.
#' @param    qc.gex            A numeric vector. Include cells and genes detected in 
#'                             at least this fraction of genes (first numeric) and 
#'                             cells (second numeric), respectively. Default is a 
#'                             minimum qc to remove zero rows and columns.
#' @param    ...               Additional arguments to be passed to \link[glasso]{glasso}. 
#'                             Otherwise, default values from the \code{glasso} 
#'                             function are used.
#' 
#' @return   a \link{genixObj} object. It can be a list of \code{genixObj}s if 
#'           \code{metadata} and \code{field} are provided.
#'
#' @details \code{constructNets} uses a glasso-based algorith to infer gene association 
#'          networks from a scRNA-seq data. 
#'
#' @seealso  \link{genixObj} for a genix object. \link[glasso]{glasso} for 
#'           details about glasso algorith.
#' 
#' @examples
#' \dontrun{
#' # getting the data:
#' data("genixDemoData")
#' 
#' # construct a graph per group:
#' rslts <- constructNets(gex=genixDemoData$gex, glasso.rho=1.0, 
#'                        metadata=genixDemoData$metadata, 
#'                        field="group", scale.gex=TRUE, 
#'                        qc.gex=c(0.01,0.5), maxit=1e3)
#' }
#' 
#' @export
constructNets <- function(gex, glasso.rho=1.0,  
                         metadata=NULL, field=NULL, 
                         scale.gex=FALSE, qc.gex=c(0,0), ...) {
  message("+++++++++++++++++++++++++++++++++++")
  start_time <- Sys.time()
  # ===================================
  # gex <- as.matrix(gex)
  # ===================================
  ## seperate data-matrices by cell type
  if (!is.null(field)) {
    message("*** ", paste("split gex by", field, "..."))
    groups <- metadata[[field]]  
  } else {
    groups <- rep("grp", ncol(gex))
  }
  verbose <- length(unique(groups)) != 1
  # ===================================
  # row -> cells
  # column -> genes
  Q <- lapply(unique(groups), function(x) {
    # subset with filed
    y <- Matrix::Matrix(gex[, which(groups == x)])
    if (verbose) {
      message("*** ", paste("A matrix of group", x, "with", dim(y)[1], "genes and", dim(y)[2], "cells received."))
    } else { 
      message("*** ", paste("A matrix with", dim(y)[1], "genes and", dim(y)[2], "cells received."))
    }
    if (scale.gex) { y <- scaleGex(y) }
    y <- qcGex(gex=y, 
               min_gene=floor(qc.gex[1]*dim(y)[1]), 
               min_cell=floor(qc.gex[2]*dim(y)[2]))$gex 
    if (verbose) {
      message("*** ", paste("A matrix of group", x, "with", dim(y)[1], "genes and", dim(y)[2], "cells proceed."))
    } else { 
      message("*** ", paste("A matrix with", dim(y)[1], "genes and", dim(y)[2], "cells proceed."))
    }
    y <- Matrix::t(y)  
    return(y)
  })
  names(Q) <- unique(groups)
  # ===================================
  ## Graphical LASSO
  # rho: the regularisation parameter that controls the degree of sparsity in 
  # the resulting inverse covariance matrix. Higher values lead to greater sparsity.
  # The covariance matrix is always both symmetric and positive semi-definite.
  # row -> genes
  # column -> genes
  glasso_res <- lapply(Q, function(x){
    if (verbose) {
      message("*** ", paste("graphical lasso", names(Q)[parent.frame()$i[]], "..."))
      message("*** ", paste0(dim(x)[2], " genes; ", dim(x)[1], " cells"))
    } else { 
      message("*** ", paste("graphical lasso ..."))  
    }
    glasso::glasso(s=cov(as.matrix(x)), rho=glasso.rho, ...)
  })
  # ===================================
  # TEST
  # glasso_res <- lapply(Q, function(x){ # [c("B", "TNK")]
  #   cat(paste(names(Q)[parent.frame()$i[]]), sep = "\n")
  #   glasso::glasso(s = cov(as.matrix(x)),  # The covariance matrix is always both symmetric and positive semi-definite.
  #                  rho = glasso.rho,
  #                  thr = 1e-4, maxit = 1e3,  
  #                  nobs = NULL, zero = NULL, 
  #                  approx = FALSE, penalize.diagonal = TRUE, 
  #                  start = c("cold","warm"), w.init = NULL, wi.init = NULL, 
  #                  trace = FALSE)
  # })
  # ===================================
  # extract precision matrix, EICS (Estimated Inverse Covariance Matrix)
  icm <- lapply(glasso_res, function(x) { x$wi })
  # ===================================
  # fix row & coulmn names
  icm <- mapply(function(x,y) {
    rownames(x) <- colnames(y)
    colnames(x) <- colnames(y)
    return(x)
  }, x = icm, y = Q, SIMPLIFY = FALSE)  
  # ===================================
  # check for symmetry in the resulting inverse covariance matrix. 
  # Assymmetry can arise due to numerical computation and rounding errors.
  icm <- lapply(icm, function(x) {
    if (!isSymmetric(x)) {
      x[lower.tri(x)] <- t(x)[lower.tri(x)]  
    }      
    return(x)
  })
  # ===================================
  # Calculate the partial correlation matrix and set the terms on the diagonal to 
  # cat(paste("calculate the partial correlation matrix ..."), sep = "\n")
  ParCorrs <- lapply(icm, function(x) {
    if (verbose) {
      message("*** ", paste("partial correlation", names(icm)[parent.frame()$i[]], "...")) 
    } else { 
      message("*** ", paste("partial correlation ..."))  
    }
    n_gns <- unique(dim(x))
    y <- matrix(0, nrow = n_gns, ncol = n_gns)
    for(i in 1:(n_gns-1)) {
      for(j in (i+1):n_gns) {
        y[i,j] <- -x[i,j]/sqrt(x[i,i]*x[j,j])
      }
    }
    y[lower.tri(y)] <- t(y)[lower.tri(y)]  
    colnames(y) <- colnames(x)
    rownames(y) <- rownames(x)
    return(y)
  })
  # print(lapply(ParCorrs, function(x) dim(x)))  
  # ===================================
  # Build the network graph
  # cat(paste("build network graph ..."), sep = "\n")
  graphs <- lapply(ParCorrs, function(x) {
    if (verbose) {
      message("*** ", paste("construct graph", names(ParCorrs)[parent.frame()$i[]], "...")) 
    } else { 
      message("*** ", paste("construct graph ..."))  
    }
    grph <- igraph::graph_from_adjacency_matrix(x, mode = "undirected", weighted = TRUE)
    return(grph)
  })
  # ===================================
  if (!verbose) {
    graphs <- graphs[[1]]
    res_obj <- new("genixObj",
                   icm_graph = graphs)
  } else {
    res_obj <- list()
    for (i in 1:length(unique(groups))) {
      genix_obj <- new("genixObj",
                       icm_graph = graphs[[i]])
      res_obj[[length(res_obj) + 1]] <- genix_obj
      names(res_obj)[length(res_obj)] <- names(graphs)[i]
    }
  }
  # ===================================
  end_time <- Sys.time()
  dt <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
  message("*** ", paste("time:", dt, "min"))
  message("+++++++++++++++++++++++++++++++++++")
  # ===================================  
  # returns
  return(res_obj)
  # ===================================  
}


#' @title Compiles a gene association network inferred by \link{constructNets}
#' @description   identifies modules and hub genes 
#'
#' @param    grph             A \link{genixObj} object. 
#' @param    degree.th        A numeric value. Genes with connectivity below \code{degree.th} 
#'                            will be removed from the graph.
#' @param    hubs.th          A numeric value. Genes that are connected to 
#'                            at least \code{hubs.th} fraction of other genes
#'                            will be considered for hub detection.
#' @param    ...              Additional arguments to be passed to 
#'                            \code{cutreeDynamicTree}. Otherwise, default values 
#'                            from the \code{cutreeDynamicTree} function are used.
#' 
#' @return   a \link{compiledObj}. It can be a list of \code{compiledObj}s, 
#'           if the \code{grph} contains a list \link{genixObj}s.
#'
#' @details \code{compileNets} processes the input graphs. First, it filters out
#'          genes with connectivity lower than \code{degree.th}. The filtered graphs 
#'          then undergo hub gene identification. A gene is considered a hub if it 
#'          exhibits an outlier level of connectivity and is connected to more than 
#'          \code{hubs.th} (fraction) of the genes in the graph. An outlier is defined as a gene 
#'          with a degree exceeding at least 1.5 interquartile ranges above the 75th 
#'          percentile of all degrees. Finally, the adjacency matrix calculated from 
#'          the graph is retrieved to identify modules using the 
#'          \code{cutreeDynamicTree} function.
#'
#' @seealso  \link{constructNets} to construct a network. 
#'           \link{compiledObj} for a \code{compileNets} output. 
#'           \code{dynamicTreeCut} package for details about 
#'            \code{cutreeDynamicTree} algorithm.
#' 
#' @examples
#' \dontrun{
#'   # first run \link{constructNets} function example and pass the \code{rslts}
#'   cmpl_rslts <- compileNets(grph=rslts, degree.th=3, hubs.th=0.01, 
#'                             minModuleSize=7, deepSplit=TRUE)
#' }
#' 
#' @export
compileNets <- function(grph, degree.th=NULL, hubs.th=NULL, ...) { 
                       # min_module_size = 10, max_tree_height = 1, deep_split = TRUE
  message("+++++++++++++++++++++++++++++++++++")
  # ======================
  num_grph <- 1
  if (is.list(grph)) {
    grph <- sapply(grph, function(x) { x@icm_graph })
    num_grph <- length(grph)
  }
  res_obj <- list()
  for (i in 1:num_grph){
    if (is.list(grph)) {
      gr <- grph[[i]]
      gr_nm <- names(grph[i])
      message("*** ", paste0(gr_nm, " graph retrieved. "))
    } else {
      gr <- grph@icm_graph
    }
    # filters
    message("*** ", paste0("A graph of ", igraph::vcount(gr), " genes and ", igraph::ecount(gr), " edges received."))
    if (!is.null(degree.th)) {
      # drop vertices with degrees less than degree.th
      rm_nodes <- which(igraph::degree(gr) < degree.th)
      gr <- igraph::delete.vertices(gr, rm_nodes)
      message("*** ", paste0(length(rm_nodes), " genes removed."))
      message("*** ", paste0("A graph of ", igraph::vcount(gr), " genes and ", igraph::ecount(gr), " edges proceed."))
    }
    # ======================      
    # find hubs
    hubs <- sort(igraph::degree(gr))
    hubs <- sort(hubs[isOutlier(hubs, tail="upper")])
    if (!is.null(hubs.th)) {
      hubs <- hubs[hubs > floor(hubs.th*(igraph::vcount(gr)-1))]
    }
    # ====================== 
    # convert graph to adjancy matrix
    g_mtx <- igraph::as_adjacency_matrix(graph=gr, type="both", attr="weight",
                                         names=TRUE, sparse=FALSE)
    stopifnot(dim(g_mtx)[1] == igraph::vcount(gr))
    stopifnot(all(rownames(g_mtx) == names(igraph::V(gr))))
    # ======================
    # module detection
    # in order to separate clusters of genes which have a strong positive
    # correlation from those with a strong negative correlation, we will first
    # shift the data from the range [-1,1] to [0,1].
    y <- 0.5 + 0.5*g_mtx # signed network (range between 0 to 1) preserves sign information
    # hierarchical clustering works by comparing the distance between objects
    # instead of the similarity.
    gene_tree <- hclust(as.dist(1-y), method = "average")
    # we will use the cuttreeDynamicTree method to break apart the hc dendrogram
    # into separate modules.
    set.seed("42")
    mdls <- dynamicTreeCut::cutreeDynamicTree(dendro=gene_tree, ...) 
    # test
    # mdls <- cutreeDynamicTree(dendro=gene_tree,
    #                              maxTreeHeight = 1,
    #                              minModuleSize = 7,
    #                              deepSplit = TRUE)
    mdls <- paste0("M", mdls)
    # ======================
    # summerize modules
    sorted_mdls <- unique(mdls)[gtools::mixedorder(unique(mdls))]
    mdls_ls <- sapply(sorted_mdls, function(x){
      logik <- which(mdls == x)
      logik <- names(igraph::V(gr))[logik]
      sort(igraph::degree(gr)[logik]) 
    }, USE.NAMES = TRUE)
    # ======================
    # add sumcorr and connectivity measures
    # connectivity identifies the strenght of a netwerk. Hubs (densly related genes)
    # will appear at the right-tail of the distribution.
    igraph::V(gr)$connectivity <- igraph::degree(gr)
    igraph::V(gr)$module <- mdls
    # ======================
    from_gene <- from_name <- from_connectivity <- from_module <- NULL
    from_module <- from_hub <- NULL
    weight <- par_corr_coeff <- NULL
    to_gene <- to_name <- to_connectivity <- to_module <- NULL
    to_module <- to_hub <- NULL
    # ======================
    g_df <- igraph::as_long_data_frame(gr) %>%
      dplyr::select(from_name, from_connectivity, from_module, 
                    weight,
                    to_name, to_connectivity, to_module) %>%
      dplyr::rename(from_gene = from_name, 
                    par_corr_coeff = weight,
                    to_gene = to_name) %>%
      dplyr::mutate(from_hub = from_gene %in% names(hubs),
                    to_hub = to_gene %in% names(hubs)) %>%
      dplyr::select(from_gene, from_connectivity, from_module, from_hub,
                    par_corr_coeff,
                    to_gene, to_connectivity, to_module, to_hub)
    # checks
    stopifnot(all(unique(g_df$from_gene[g_df$from_hub]) %in% names(hubs)))
    stopifnot(all(unique(g_df$to_gene[g_df$to_hub]) %in% names(hubs)))
    stopifnot(all(names(hubs) %in% unique(c(g_df$from_gene[g_df$from_hub], g_df$to_gene[g_df$to_hub]))))
    # ======================
    # print summary
    n_modules <- length(unique(mdls)) - 1
    n_hubs <- length(hubs)
    message("*** ", paste(n_modules, "module(s) and", n_hubs, "hub(s) detected."))
    # ======================
    compiled_obj <- new("compiledObj",
                        graph = gr,
                        graph_df = g_df,
                        hubs = hubs,
                        modules = mdls_ls)
    # ======================
    if (is.list(grph)) {
      res_obj[[length(res_obj) + 1]] <- compiled_obj
      names(res_obj)[length(res_obj)] <- gr_nm
    } else {
      res_obj <- compiled_obj
    } 
  }
  message("+++++++++++++++++++++++++++++++++++")
  # ======================
  # returns
  return(res_obj)
  # ===================================  
}


#' @title Compares topological characteristics of two graphs
#' @description   identifies the variations in topological features of two graphs
#'
#' @param    grph.1           A \code{compiledObj} object. 
#' @param    grph.2           A \code{compiledObj} object. 
#' @param    n.perm           An integer. Number of times to perform permutation
#'                            test between two graphs.
#' 
#' @return   a \link{comparedObj}.
#'
#' @details \code{compareNets} calculates the normalized values of degree,
#'   closeness, and betweenness for each gene in two graphs, grph.1 and grph.2. 
#'   It then substrates the calculated topological characteristics to report deltas.
#'   Positive values indicate that the feature is larger in the grph.2. 
#'   
#'   The function also calculates the similarity between grph.1 and grph.2 using the
#'   Jaccard index, considering both shared genes and edges. It then performs an
#'   iterative gene removal process in both graphs. With each gene removal, 
#'   the similarity index between graphs is recalculated. The gene removal impact 
#'   for each gene is then calculated as the relative difference in post-removal
#'   graph similarity compared to the pre-removal.
#'   
#'   A permutation test is also developed to examine the statistical distinctiveness
#'   of the \code{grph.1} and \code{grph.2} First, the two graphs are combined (union) into a
#'   pull. Next, \code{grph.1} and \code{grph.2} are reconstructed by randomly sampling from
#'   the pull, while preserving the original edge counts and connections, with
#'   no-replacement. A distribution of jaccard similarities is generated
#'   (n=\code{n.perm}). This distribution is utilized to calculate a two-sided
#'   t-test p-value, reflecting the proportion of permuted statistics that
#'   differ from the observed jaccard index between \code{grph.1} and \code{grph.2}
#'
#' @seealso  \link{constructNets} to construct a network.
#'           \link{compiledObj} for a \link{compileNets} output. 
#' 
#' @examples
#' \dontrun{
#'   # first run \link{compileNets} function example and pass the \code{cmpl_rslts}
#'   cmpr_rslts <- compareNets(grph.1=cmpl_rslts$grp1, grph.2=cmpl_rslts$grp2, 
#'                             n.perm=25)
#' }
#' 
#' @export
compareNets <- function(grph.1, grph.2, n.perm=1000) {
  message("+++++++++++++++++++++++++++++++++++")
  # ===================================
  grph1 <- grph.1@graph 
  hubs1 <- names(grph.1@hubs)
  
  grph2 <- grph.2@graph
  hubs2 <- names(grph.2@hubs)
  # ===================================
  message("*** grph.1: ", paste("A graph of", igraph::vcount(grph1), "genes and", 
                                 igraph::ecount(grph1), "edges received."))
  message("*** grph.2: ", paste("A graph of", igraph::vcount(grph2), "genes and", 
                                 igraph::ecount(grph2), "edges received."))
  # ===================================
  g.intersect <- igraph::intersection(grph1, grph2, 
                                      byname=TRUE, keep.all.vertices=FALSE)
  # message(paste("*** intersection grph.1 and grph.2:", 
  #               "A graph of", igraph::vcount(g.intersect), 
  #               "genes and", igraph::ecount(g.intersect), "edges."))
  # ===================================
  g.union <- igraph::union(grph1, grph2, byname=TRUE)
  # message(paste("*** union grph.1 and grph.2:", 
  #               "A graph of", igraph::vcount(g.union), 
  #               "genes and", igraph::ecount(g.union), "edges."))
  # ===================================
  sim_index <- igraph::ecount(g.intersect)/igraph::ecount(g.union)
  # ===================================
  # approximate permutation test
  message(paste("*** graphs distinctiveness test..."))
  ecount_ga <- igraph::ecount(grph1)
  ecount_gb <- igraph::ecount(grph2)
  prmts <- replicate(n.perm, {
    ga_shfl <- igraph::subgraph.edges(graph=g.union, 
                              eids=sample(igraph::E(g.union), ecount_ga, replace=FALSE),
                              delete.vertices=TRUE)
    gb_shfl <- igraph::subgraph.edges(graph=g.union, 
                              eids=sample(igraph::E(g.union), ecount_gb, replace=FALSE),
                              delete.vertices=TRUE)
    g.intersect_shfl <- igraph::intersection(ga_shfl, gb_shfl, 
                                             byname=TRUE, keep.all.vertices=FALSE)
    g.union_shfl <- igraph::union(ga_shfl, gb_shfl, byname=TRUE)
    igraph::ecount(g.intersect_shfl)/igraph::ecount(g.union_shfl)
  })
  # null hypothesis: two graphs are similar (no difference between conditions)
  pvalue <- t.test(x = prmts, mu = sim_index, alternative = "two.sided")$p.value
  # ===================================
  a_genes <- igraph::V(grph1)$name
  b_genes <- igraph::V(grph2)$name
  # ===================================
  conf_mtx <- twoSetsCommp(A=a_genes, B=b_genes)
  in_A <- setdiff(a_genes, b_genes)
  in_B <- setdiff(b_genes, a_genes)
  in_AB <- intersect(a_genes, b_genes)
  in_all <- union(a_genes, b_genes)
  # ===================================
  message("*** gene removal impact...")
  reverse_comp <- sapply(in_all, function(gn){
    if (gn %in% a_genes) { ga <- igraph::delete.vertices(grph1, gn) } else { ga <- grph1 }
    if (gn %in% b_genes) { gb <- igraph::delete.vertices(grph2, gn) } else { gb <- grph2 }
    g.intersect <- igraph::intersection(ga, gb, byname=TRUE, keep.all.vertices=FALSE)
    g.union <- igraph::union(ga, gb, byname=TRUE)
    return(igraph::ecount(g.intersect)/igraph::ecount(g.union))
  }) 
  igraph::E(grph1)$weight <- rep(1, igraph::ecount(grph1))
  igraph::E(grph2)$weight <- rep(1, igraph::ecount(grph2))
  # ===================================
  gene <- jaccard <- delta_jaccard <- NULL
  degree.1 <- closeness.1 <- betweenness.1 <- NULL
  degree.2 <- closeness.2 <- betweenness.2 <- NULL
  # ===================================
  dfp2 <- data.frame(gene = names(reverse_comp),
                     jaccard = as.vector(reverse_comp)) %>%
    dplyr::mutate(degree.1 = igraph::degree(grph1, normalized=TRUE)[gene],
                  closeness.1 = igraph::closeness(grph1, normalized=TRUE, weights=abs(igraph::E(grph1)$weight))[gene],
                  betweenness.1 = igraph::betweenness(grph1, normalized=TRUE, weights=abs(igraph::E(grph1)$weight))[gene],
                  degree.2 = igraph::degree(grph2, normalized=TRUE)[gene],
                  closeness.2 = igraph::closeness(grph2, normalized=TRUE, weights=abs(igraph::E(grph2)$weight))[gene],
                  betweenness.2 = igraph::betweenness(grph2, normalized=TRUE, weights=abs(igraph::E(grph2)$weight))[gene],
                  delta_jaccard = (jaccard - sim_index)/sim_index)
  dfp2[is.na(dfp2)] <- 0
  dfp2$delta_degree <- dfp2$degree.2 - dfp2$degree.1
  dfp2$delta_closeness <- dfp2$closeness.2 - dfp2$closeness.1
  dfp2$delta_betweenness <- dfp2$betweenness.2 - dfp2$betweenness.1
  # ===================================
  if (!is.null(hubs1) & !is.null(hubs2)) {
    dfp2$type <- "common"
    dfp2$type[dfp2$gene %in% hubs1] <- "hub.1"
    dfp2$type[dfp2$gene %in% hubs2] <- "hub.2"
    dfp2$type[dfp2$gene %in% intersect(hubs1, hubs2)] <- "hub.12"
  }
  # ===================================
  message("+++++++++++++++++++++++++++++++++++")
  # returns
  compared_obj <- new("comparedObj",
                      ji = sim_index,
                      pvalue = pvalue,
                      permutations = prmts,
                      deltas = dfp2,
                      inA = in_A,
                      inB = in_B,
                      inAB = in_AB,
                      inAB_mtx = conf_mtx)
  return(compared_obj)
  # ===================================
}


#' @title Topological variant gene identification
#' @description   filter genes with topological feastures above a threshold.
#'
#' @param    obj       a \link{comparedObj} object.
#' @param    feature   a character. The feature by which to filter the data.
#' @param    sigma     a numeric. The standard deviation coefficient by which to 
#'                     filter the data.
#'   
#' @return   A data.frame containing gene names that have passed the filter of
#'   \code{sigma}*\code{sd}(\code{feature}). The output data.frame is sorted by
#'   the values of \code{feature} in descending order.
#'
#' @examples
#' \dontrun{
#' # first run \link{compareNets} function example and pass the \code{cmpr_rslts} 
#' 
#' x <- fetchTVGs(cmpr_rslts, feature="delta_degree", sigma=1)
#'}
#' 
#' @export
fetchTVGs <- function(obj, feature, sigma=1) {
  # Hack for visibility of dplyr variables
  . <- NULL
  x <- NULL
  
  thrshld <- sigma*sd(obj@deltas[[feature]])
  message(paste("*** fitering genes with absolute value larger than", round(thrshld, 4), "in", feature, "feature."))
  x <- obj@deltas %>%
    dplyr::filter(abs(!!as.symbol(feature)) > thrshld) %>%
    dplyr::arrange(desc(!!as.symbol(feature)))
  return(x)
}


#### Summary functions ####


#' @title plots inferred network characteristics
#' @description   Displays the distribution of degrees in a graph.
#'
#' @param    obj       a \link{compiledObj} object.
#'   
#' @return   A ggplot object.
#'
#' @examples
#' \dontrun{
#' # first run  \link{compileNets} function example and pass the \link{cmpl_rslts} 
#' 
#' # plot degree distribution:
#' plotNetDegree(cmpl_rslts$grp2)
#'}
#' 
#' @export
plotNetDegree <- function(obj) {
  # Hack for visibility of dplyr variables
  . <- NULL
  Freq <- x <- .x <- y <- NULL
  # ======================
  # Transitivity is the overall probability for the network to have adjacent 
  # nodes interconnected, thus revealing the existence of tightly connected 
  # communities (or clusters, subgroups, cliques). The transitivity of a 
  # fully-connected graph is equal to 1.
  # type: The type of the transitivity to calculate. Possible values:
  # "global": The global transitivity of an undirected graph. This is simply the 
  # ratio of the count of triangles and connected triples in the graph. 
  g.degrees <- igraph::degree(obj@graph)
  g.degree.hist <- data.frame(g.degrees = as.numeric(names(table(g.degrees))),
                              Freq =as.vector(table(g.degrees)))
  p <- ggplot(g.degree.hist, 
              aes(x = g.degrees, y = Freq)) +
    theme_bw() + 
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 15, face = "bold")) +
    ggtitle("Degree Distribution", subtitle = paste0("Graph transitivity: ", round(igraph::transitivity(obj@graph), 2))) +
    xlab("Connectivity") + ylab("Frequency") + 
    geom_point(size = 2, shape = 1) +
    geom_smooth(size = 0.35, method = 'loess', formula = y ~ x) +
    scale_x_continuous(trans = scales::log10_trans(),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    scale_y_continuous(trans = scales::log10_trans(),
                       breaks = scales::trans_breaks("log10", function(x) 10^x),
                       labels = scales::trans_format("log10", scales::math_format(10^.x)))
  return(p)
}


#' @title graphs distinctiveness
#' @description  Displays the distribution of Jaccard index similarities 
#'
#' @param    obj       a \link{comparedObj} object.
#'   
#' @return   A ggplot object.
#'
#' @seealso  \code{compareNets} to generate \link{comparedObj}. 
#'           
#' @examples
#' \dontrun{
#' # first run \link{compareNets} function example and pass the \code{cmpr_rslts} 
#' 
#' # plot permutation distribution:
#'  plotPermutes(cmpr_rslts)
#'}
#' 
plotPermutes <- function(obj){
  # Hack for visibility of dplyr variables
  . <- NULL
  x <- NULL
  # ======================
  stopifnot(class(obj)[1] == "comparedObj")
  df <- data.frame(x = obj@permutations)
  p <- ggplot(df, aes(x=x)) +
    theme_bw() + 
    theme(axis.text = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 15, face = "bold")) +
    ggtitle(paste("p-value:", round(obj@pvalue, 3))) +
    xlab("Similarity") + ylab("Count") +
    geom_histogram(bins = 100, fill="steelblue", color="white") +
    geom_vline(xintercept = obj@ji, color = "firebrick", linetype = 2)
  return(p)
}


#### Helper functions ####


# Helper function for gex qc
qcGex <- function(gex, min_gene=0, min_cell=0) {
  dims_i <- dim(gex)
  while (any(Matrix::colSums(gex != 0) < min_gene) | any(Matrix::rowSums(gex != 0) < min_cell)) {
    # filtered out cells for which fewer than min_gene genes were detected,
    # and genes that were expressed in fewer than min_cell cells.
    gex <- gex[, Matrix::colSums(gex != 0) >= min_gene ]
    gex <- gex[Matrix::rowSums(gex != 0) >= min_cell , ]
  }
  dims_f <- dim(gex)
  # ===================================
  message("*** ", paste(dims_i[1] - dims_f[1], "gene(s) and", dims_i[2] - dims_f[2], "cell(s) removed."))
  # ===================================
  # returns
  return(list("gex" = gex))
}


# scale input gene expression matrix
scaleGex <- function (gex) {
  xx = NULL
  if (!is.null(colnames(gex))) 
    xx <- colnames(gex)
  if (class(gex)[1] == "matrix") {
    m <- Matrix::Matrix(0, ncol(gex), ncol(gex))
    tots_use <- Matrix::colSums(gex)
    target_mean <- mean(tots_use)
    diag(m) <- target_mean/tots_use
    gex <- as.matrix(gex %*% m)
    if (!is.null(xx)) {
      colnames(gex) <- xx
    }
  }
  else {
    m <- Matrix::Matrix(0, ncol(gex), ncol(gex))
    tots_use <- Matrix::colSums(gex)
    target_mean <- mean(tots_use)
    diag(m) <- target_mean/tots_use
    gex <- Matrix::Matrix(gex %*% m, sparse = TRUE)
    if (!is.null(xx)) {
      colnames(gex) <- xx
    }
  }
  return(gex)
}


# Helper function for identifying outliers
isOutlier <- function(x, tail = "both") {
  if (tail == "both") {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))  
  } else if (tail == "lower") {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x))  
  } else if (tail == "upper") {
    return(x > quantile(x, 0.75) + 1.5 * IQR(x))  
  }
}


# comapre genes in two sets
twoSetsCommp <- function(A,B){
  both    <-  union(A,B)
  inA     <-  both %in% A
  inB     <-  both %in% B
  return(table(inA,inB))
}
