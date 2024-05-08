#' @rdname getLineages
#'
#' @description This function constructs the minimum spanning tree(s) on
#'   clusters of cells, the first step in Slingshot's trajectory inference
#'   procedure. Paths through the MST from an origin cluster to leaf node
#'   clusters are interpreted as lineages.
#'
#' @param data a data object containing the matrix of coordinates to be used for
#'   lineage inference. Supported types include \code{matrix},
#'   \code{\link{SingleCellExperiment}}, \code{\link{SlingshotDataSet}}, and
#'   \code{\link[TrajectoryUtils]{PseudotimeOrdering}}.
#' @param clusterLabels each cell's cluster assignment. This can be a single
#'   vector of labels, or a \code{#cells} by \code{#clusters} matrix
#'   representing weighted cluster assignment. Either representation may
#'   optionally include a \code{"-1"} group meaning "unclustered."
#' @param reducedDim (optional) the dimensionality reduction to be used. Can be
#'   a matrix or a character identifying which element of
#'   \code{reducedDim(data)} is to be used. If multiple dimensionality
#'   reductions are present and this argument is not provided, the first element
#'   will be used by default.
#' @param start.clus (optional) character, indicates the starting cluster(s)
#'   from which lineages will be drawn.
#' @param end.clus (optional) character, indicates which cluster(s) will be
#'   forced to be leaf nodes in the graph.
#' @param dist.method (optional) character, specifies the method for calculating
#'   distances between clusters. Default is \code{"slingshot"}, see
#'   \code{\link[TrajectoryUtils]{createClusterMST}} for details.
#' @param use.median logical, whether to use the median (instead of mean) when
#'   calculating cluster centroid coordinates.
#' @param omega (optional) numeric or logical, this granularity parameter
#'   determines the distance between every real cluster and the artificial
#'   cluster, \code{.OMEGA}. In practice, this makes \code{omega} the maximum
#'   allowable distance between two connected clusters. By default, \code{omega
#'   = Inf}. If \code{omega = TRUE}, the maximum edge length will be set to the
#'   median edge length of the unsupervised MST times a scaling factor
#'   (\code{omega_scale}, default \code{= 1.5}). This value is provided as a
#'   potentially useful rule of thumb for datasets with outlying clusters or
#'   multiple, distinct trajectories. See \code{outgroup} in
#'   \code{\link[TrajectoryUtils]{createClusterMST}}.
#' @param omega_scale (optional) numeric, scaling factor to use when \code{omega
#'   = TRUE}. The maximum edge length will be set to the median edge length of
#'   the unsupervised MST times \code{omega_scale} (default \code{= 3}). See
#'   \code{outscale} in \code{\link[TrajectoryUtils]{createClusterMST}}.
#' @param times numeric, vector of external times associated with either
#'   clusters or cells. See \code{\link[TrajectoryUtils]{defineMSTPaths}} for
#'   details.
#'
#' @details Given a reduced-dimension data matrix \code{n} by \code{p} and a set
#'   of cluster identities (potentially including a \code{"-1"} group for
#'   "unclustered"), this function infers a tree (or forest) structure on the
#'   clusters. This work is now mostly handled by the more general function,
#'   \code{\link[TrajectoryUtils]{createClusterMST}}.
#'
#' @details The graph of this structure is learned by fitting a (possibly
#'   constrained) minimum-spanning tree on the clusters, plus the artificial
#'   cluster, \code{.OMEGA}, which is a fixed distance away from every real
#'   cluster. This effectively limits the maximum branch length in the MST to
#'   the chosen distance, meaning that the output may contain multiple trees.
#'
#' @details Once the graph is known, lineages are identified in
#'   any tree with at least two clusters. For a given tree, if there is an
#'   annotated starting cluster, every possible path out of a starting cluster
#'   and ending in a leaf that isn't another starting cluster will be returned.
#'   If no starting cluster is annotated, one will be chosen by a heuristic
#'   method, but this is not recommended.
#'
#' @return An object of class \code{\link{PseudotimeOrdering}}. Although the
#'   final pseudotimes have not yet been calculated, the assay slot of this
#'   object contains two elements: \code{pseudotime}, a matrix of \code{NA}
#'   values; and \code{weights}, a preliminary matrix of lineage assignment
#'   weights. The \code{reducedDim} and \code{clusterLabels} matrices will be
#'   stored in the \code{\link[TrajectoryUtils]{cellData}}. Additionally, the
#'   \code{metadata} slot will contain an \code{\link[igraph]{igraph}} object
#'   named \code{mst}, a list of parameter values named \code{slingParams}, and
#'   a list of lineages (ordered sets of clusters) named \code{lineages}.
#'
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- getLineages(rd, cl, start.clus = '1')
#'
#' # plotting
#' sds <- as.SlingshotDataSet(pto)
#' plot(rd, col = cl, asp = 1)
#' lines(sds, type = 'l', lwd = 3)
#'
#' @export
#'
#' @import TrajectoryUtils
#'
setMethod(f = "getLineages",
    signature = signature(data = "matrix",
        clusterLabels = "matrix"),
    definition = function(data, clusterLabels, reducedDim = NULL,
        start.clus = NULL, end.clus = NULL,
        dist.method = "slingshot", use.median = FALSE,
        omega = FALSE, omega_scale = 1.5,
        times = NULL,
        distmat = NULL, ...){

        X <- as.matrix(data)
        clusterLabels <- as.matrix(clusterLabels)
        ####################
        ### CHECKS
        ####################
        if(nrow(X)==0){
            stop('reducedDim has zero rows.')
        }
        if(ncol(X)==0){
            stop('reducedDim has zero columns.')
        }
        if(nrow(X) != nrow(clusterLabels)){
            stop('nrow(data) must equal nrow(clusterLabels).')
        }
        if(any(is.na(X))){
            stop('reducedDim cannot contain missing values.')
        }
        if(!all(apply(X,2,is.numeric))){
            stop('reducedDim must only contain numeric values.')
        }
        if (is.null(rownames(X)) &
                is.null(rownames(clusterLabels))) {
            rownames(X) <- paste('Cell', seq_len(nrow(X)), sep = '-')
            rownames(clusterLabels) <-
                paste('Cell', seq_len(nrow(X)), sep = '-')
        }
        if(is.null(colnames(X))){
            colnames(X) <- paste('Dim',seq_len(ncol(X)),sep='-')
        }
        if(is.null(colnames(clusterLabels))) {
            colnames(clusterLabels) <- seq_len(ncol(clusterLabels))
        }
        if(any(colnames(clusterLabels) == "")){
            colnames(clusterLabels)[colnames(clusterLabels)==""] <-
                which(colnames(clusterLabels)=="")
        }
        if(any(rownames(X)=='')){
            miss.ind <- which(rownames(X) == '')
            rownames(X)[miss.ind] <- paste('Cell',miss.ind,sep='-')
        }
        if(any(colnames(X)=='')){
            miss.ind <- which(colnames(X) == '')
            colnames(X)[miss.ind] <- paste('Dim',miss.ind,sep='-')
        }
        if(is.null(rownames(clusterLabels)) &
                !is.null(rownames(X))){
            rownames(clusterLabels) <- rownames(X)
        }
        if(is.null(rownames(X)) &
                !is.null(rownames(clusterLabels))){
            rownames(X) <- rownames(clusterLabels)
        }
        if(any(rowSums(clusterLabels)>1)){
            rs <- rowSums(clusterLabels)
            clusterLabels <- clusterLabels / rs
        }
        if(any(colSums(clusterLabels)==0)){
            clusterLabels <- clusterLabels[, colSums(clusterLabels)!=0,
                drop = FALSE]
        }
        if(length(omega) > 1){
            stop('omega must have length 1')
        }
        if(is.na(as.numeric(omega))){
            stop('omega must be logical or numeric')
        }
        if(omega < 0){
            stop('omega must be non-negative')
        }

        # set up, remove unclustered cells (-1's)
        clusterLabels <- clusterLabels[, colnames(clusterLabels) != -1,
            drop = FALSE]
        clusters <- colnames(clusterLabels)
        nclus <- length(clusters)
        if(!is.null(start.clus)){
            start.clus <- as.character(start.clus)
        }
        if(!is.null(end.clus)){
            end.clus <- as.character(end.clus)
        }

        ### make the MST / forest (multiple MSTs)
        if(nclus == 1){
            dmat <- matrix(0)
            rownames(dmat) <- colnames(dmat) <- clusters
            g <- igraph::graph_from_adjacency_matrix(dmat,
                                                     mode = "undirected",
                                                     weighted = TRUE)
            igraph::V(g)$coordinates <- list(clusters = colMeans(X))

            lineages <- list('Lineage1' = clusters)
        }else{
            use <- which(rowSums(clusterLabels) > 0)
            g <- createClusterMST_adapt(x = X[use, ,drop=FALSE],
                                  clusters = clusterLabels[use, ,drop=FALSE],
                                  outgroup = omega, outscale = omega_scale,
                                  endpoints = end.clus,
                                  dist.method = dist.method,
                                  use.median = use.median,
                                  distmat = distmat, ...)


            # select root nodes (one per connected component of g)
            forest <- igraph::decompose(g)
            starts <- vapply(forest, function(tree){
                if(length(igraph::V(tree)) == 1){
                    return(names(igraph::V(tree)))
                }
                if(any(start.clus %in% names(igraph::V(tree)))){
                    return(start.clus[start.clus %in%
                                          names(igraph::V(tree))][1])
                }
                # otherwise, pick root based on highest average length of
                # the resulting lineages (~parsimony, maximizing shared parts)
                adj <- igraph::as_adjacency_matrix(tree, sparse = FALSE)
                leaves <- rownames(adj)[rowSums(adj) == 1]
                avg.lineage.length <- vapply(leaves,function(l){
                    ends <- leaves[leaves != l]
                    paths <- igraph::shortest_paths(tree, from = l, to = ends,
                                                    mode = 'out',
                                                    output = 'vpath')$vpath
                    mean(vapply(paths, length, 0))
                }, 0)
                return(names(avg.lineage.length)[which.max(avg.lineage.length)])
            }, FUN.VALUE = '')

            lineages <- TrajectoryUtils::defineMSTPaths(g, roots = starts,
                                                        times = times,
                                                    clusters = clusterLabels,
                                                    use.median = use.median)

            # sort by number of clusters included
            lineages <- lineages[order(vapply(lineages, length, 0),
                                       decreasing = TRUE)]
            names(lineages) <- paste('Lineage',seq_along(lineages),sep='')
        }

        lineageControl <- list()
        first <- unique(vapply(lineages,function(l){ l[1] },''))
        last <- unique(vapply(lineages,function(l){ l[length(l)] },''))

        lineageControl$start.clus <- first
        lineageControl$end.clus <- last

        start.given <- first %in% start.clus
        end.given <- last %in% end.clus
        lineageControl$start.given <- start.given
        lineageControl$end.given <- end.given

        lineageControl$omega <- omega
        lineageControl$omega_scale <- omega_scale

        # cells x lineages weights matrix
        W <- vapply(seq_along(lineages),function(l){
            rowSums(clusterLabels[, lineages[[l]], drop = FALSE])
        }, rep(0,nrow(X))) # weighting matrix
        rownames(W) <- rownames(X)
        colnames(W) <- names(lineages)

        # empty pseudotime matrix (to be filled by getCurves)
        pst <- W
        pst[,] <- NA

        out <- PseudotimeOrdering(pathStats = list(pseudotime = pst,
                                                   weights = W),
                                  metadata = list(lineages = lineages,
                                                  mst = g,
                                                  slingParams = lineageControl))
        cellData(out)$reducedDim <- X
        cellData(out)$clusterLabels <- clusterLabels

        validObject(out)
        return(out)
    }
)

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", clusterLabels = "character"),
    definition = function(data, clusterLabels, ...){
        # CHECKS
        clusterLabels <- as.character(clusterLabels)
        X <- as.matrix(data)
        if(nrow(X) != length(clusterLabels)){
            stop('nrow(data) must equal length(clusterLabels).')
        }
        if(any(is.na(clusterLabels))){
            message("Cluster labels of 'NA' being treated as unclustered.")
            clusterLabels[is.na(clusterLabels)] <- '-1'
        }

        # convert clusterLabels into cluster weights matrix
        clusters <- unique(clusterLabels)
        clusWeight <- vapply(clusters,function(clID){
            as.numeric(clusterLabels == clID)
        },rep(0,nrow(X)))
        colnames(clusWeight) <- clusters
        return(getLineages(data = data, clusterLabels = clusWeight, ...))
    }
)

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix", clusterLabels = "ANY"),
    definition = function(data, clusterLabels, ...){
        if(missing(clusterLabels)){
            message('No cluster labels provided. Continuing with ',
                'one cluster.')
            clusterLabels <- rep('1', nrow(data))
        }
        if(! any(c(length(clusterLabels), nrow(clusterLabels)) ==
                nrow(data))){
            stop("clusterLabels must have length or number of rows equal',
                'to nrow(data).")
        }
        return(getLineages(data = data, clusterLabels = clusterLabels, ...))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "SlingshotDataSet",
        clusterLabels = "ANY"),
    definition = function(data, clusterLabels, ...){
        return(getLineages(data = reducedDim(data),
            clusterLabels = slingClusterLabels(data), ...))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
          signature = signature(data = "PseudotimeOrdering",
                                clusterLabels = "ANY"),
          definition = function(data, clusterLabels, ...){
              return(getLineages(data = cellData(data)$reducedDim,
                                 clusterLabels = cellData(data)$clusterLabels,
                                 ...))
          })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "data.frame",
        clusterLabels = "ANY"),
    definition = function(data, clusterLabels, ...){
        RD <- as.matrix(data)
        rownames(RD) <- rownames(data)
        return(getLineages(data = RD, clusterLabels = clusterLabels, ...))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix",
        clusterLabels = "numeric"),
    definition = function(data, clusterLabels, ...){
        return(getLineages(data = data,
            clusterLabels = as.character(clusterLabels), ...))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
    signature = signature(data = "matrix",
        clusterLabels = "factor"),
    definition = function(data, clusterLabels, ...){
        return(getLineages(data = data,
            clusterLabels = as.character(clusterLabels), ...))
    })

#' @rdname getLineages
#' @export
setMethod(f = "getLineages",
          signature = signature(data = "SingleCellExperiment"),
          definition = function(data, clusterLabels, reducedDim = NULL, ...){
            # SETUP
            # determine the cluster labels and reducedDim matrix
            if(is.null(reducedDim)){
              if(length(reducedDims(data))==0){
                stop('No dimensionality reduction found.')
              }else{
                message('Dimensionality reduction not explicitly ',
                        'chosen. Continuing with ',
                        names(reducedDims(data))[1])
                rd <- reducedDims(data)[[1]]
              }
            }
            if(length(reducedDim)==1){
              if(reducedDim %in% names(reducedDims(data))){
                rd <- reducedDims(data)[[as.character(reducedDim)]]
              }else{
                stop(reducedDim,' not found in reducedDims(data).')
              }
            }else{
              if(!is.null(dim(reducedDim))){
                rd <- reducedDim
                reducedDims(data)$slingReducedDim <- reducedDim
              }
            }

            if(missing(clusterLabels)){
              message('No cluster labels provided. Continuing with one ',
                      'cluster.')
              cl <- rep('1', nrow(rd))
            }else{
              if(length(clusterLabels)==1){
                if(clusterLabels %in% colnames(colData(data))){
                  cl <- colData(data)[[as.character(clusterLabels)]]
                }else{
                  stop(clusterLabels,' not found in colData(data).')
                }
              }
              if(length(clusterLabels)>1){
                if(!is.null(dim(clusterLabels)) &&
                   length(dim(clusterLabels)) > 1 &&
                   all(dim(clusterLabels) > 1)){
                  cl <- as.matrix(clusterLabels)
                  colnames(cl) <- paste0('sling_c',seq_len(ncol(cl)))
                }else{
                  cl <- as.character(clusterLabels)
                }
              }
            }
            # run slingshot
            pto <- getLineages(data = rd, clusterLabels = cl,
                              reducedDim = NULL, ...)
            # combine getLineages output with SCE
            sce <- data
            colData(sce)$slingshot <- pto
            return(sce)
          })


