#' @rdname slingshot
#'
#' @description Perform trajectory inference by (1) identifying lineage
#'   structure with a cluster-based minimum spanning tree, and (2) constructing
#'   smooth representations of each lineage using simultaneous principal curves.
#'   This function wraps the \code{\link{getLineages}} and
#'   \code{\link{getCurves}} functions and is the primary function of the
#'   \code{slingshot} package.
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
#' @param omega (optional) numeric, this granularity parameter determines the
#'   distance between every real cluster and the artificial cluster,
#'   \code{.OMEGA}. In practice, this makes \code{omega} the maximum allowable
#'   distance between two connected clusters. By default, \code{omega = Inf}. If
#'   \code{omega = TRUE}, the maximum edge length will be set to the median edge
#'   length of the unsupervised MST times a scaling factor (\code{omega_scale},
#'   default \code{= 3}). This value is provided as a potentially useful rule of
#'   thumb for datasets with outlying clusters or multiple, distinct
#'   trajectories. See \code{outgroup} in
#'   \code{\link[TrajectoryUtils]{createClusterMST}}.
#' @param omega_scale (optional) numeric, scaling factor to use when \code{omega
#'   = TRUE}. The maximum edge length will be set to the median edge length of
#'   the unsupervised MST times \code{omega_scale} (default \code{= 1.5}). See
#'   \code{outscale} in \code{\link[TrajectoryUtils]{createClusterMST}}.
#' @param times numeric, vector of external times associated with either
#'   clusters or cells. See \code{\link[TrajectoryUtils]{defineMSTPaths}} for
#'   details.
#' @param shrink logical or numeric between 0 and 1, determines whether and how
#'   much to shrink branching lineages toward their average prior to the split
#'   (default \code{= TRUE}).
#' @param extend character, how to handle root and leaf clusters of lineages
#'   when constructing the initial, piece-wise linear curve. Accepted values are
#'   \code{'y'} (default), \code{'n'}, and \code{'pc1'}. See 'Details' for more.
#' @param reweight logical, whether to allow cells shared between lineages to be
#'   reweighted during curve fitting. If \code{TRUE} (default), cells shared
#'   between lineages will be iteratively reweighted based on the quantiles of
#'   their projection distances to each curve. See 'Details' for more.
#' @param reassign logical, whether to reassign cells to lineages at each
#'   iteration. If \code{TRUE} (default), cells will be added to a lineage when
#'   their projection distance to the curve is less than the median distance for
#'   all cells currently assigned to the lineage. Additionally, shared cells
#'   will be removed from a lineage if their projection distance to the curve is
#'   above the 90th percentile and their weight along the curve is less than
#'   \code{0.1}.
#' @param thresh numeric, determines the convergence criterion. Percent change
#'   in the total distance from cells to their projections along curves must be
#'   less than \code{thresh}. Default is \code{0.001}, similar to
#'   \code{\link[princurve]{principal_curve}}.
#' @param maxit numeric, maximum number of iterations (default \code{= 15}), see
#'   \code{\link[princurve]{principal_curve}}.
#' @param stretch numeric factor by which curves can be extrapolated beyond
#'   endpoints. Default is \code{2}, see
#'   \code{\link[princurve]{principal_curve}}.
#' @param approx_points numeric, whether curves should be approximated by a
#'   fixed number of points. If \code{FALSE} (or 0), no approximation will be
#'   performed and curves will contain as many points as the input data. If
#'   numeric, curves will be approximated by this number of points (default
#'   \code{= 150} or \code{#cells}, whichever is smaller). See 'Details' and
#'   \code{\link[princurve]{principal_curve}} for more.
#' @param smoother choice of scatter plot smoother. Same as
#'   \code{\link[princurve]{principal_curve}}, but \code{"lowess"} option is
#'   replaced with \code{"loess"} for additional flexibility.
#' @param shrink.method character denoting how to determine the appropriate
#'   amount of shrinkage for a branching lineage. Accepted values are the same
#'   as for \code{kernel} in \code{\link{density}} (default is \code{"cosine"}),
#'   as well as \code{"tricube"} and \code{"density"}. See 'Details' for more.
#' @param allow.breaks logical, determines whether curves that branch very close
#'   to the origin should be allowed to have different starting points.
#' @param ... Additional parameters to pass to scatter plot smoothing function,
#'   \code{smoother}.
#'
#'
#' @details Given a reduced-dimensional data matrix \code{n} by \code{p} and a
#'   vector of cluster labels (or matrix of soft cluster assignments,
#'   potentially including a \code{-1} label for "unclustered"), this function
#'   performs trajectory inference using a cluster-based minimum spanning tree
#'   on the clusters and simultaneous principal curves for smooth, branching
#'   paths.
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
#' @details When there is only a single lineage, the curve-fitting algorithm is
#'   nearly identical to that of \code{\link[princurve]{principal_curve}}. When
#'   there are multiple lineages and \code{shrink > 0}, an additional step
#'   is added to the iterative procedure, forcing curves to be similar in the
#'   neighborhood of shared points (ie., before they branch).
#'
#' @details The \code{approx_points} argument, which sets the number of points
#'   to be used for each curve, can have a large effect on computation time. Due
#'   to this consideration, we set the default value to \code{150} whenever the
#'   input dataset contains more than that many cells. This setting should help
#'   with exploratory analysis while having little to no impact on the final
#'   curves. To disable this behavior and construct curves with the maximum
#'   number of points, set \code{approx_points = FALSE}.
#'
#' @details The \code{extend} argument determines how to construct the
#'   piece-wise linear curve used to initiate the recursive algorithm. The
#'   initial curve is always based on the lines between cluster centers and if
#'   \code{extend = 'n'}, this curve will terminate at the center of the
#'   endpoint clusters. Setting \code{extend = 'y'} will allow the first and
#'   last segments to extend beyond the cluster center to the orthogonal
#'   projection of the furthest point. Setting \code{extend = 'pc1'} is similar
#'   to \code{'y'}, but uses the first principal component of the cluster to
#'   determine the direction of the curve beyond the cluster center. These
#'   options typically have limited impact on the final curve, but can
#'   occasionally help with stability issues.
#'
#' @details When \code{shink = TRUE}, we compute a percent shrinkage curve,
#'   \eqn{w_l(t)}, for each lineage, a non-increasing function of pseudotime
#'   that determines how much that lineage should be shrunk toward a shared
#'   average curve. We set \eqn{w_l(0) = 1} (complete shrinkage), so that the
#'   curves will always perfectly overlap the average curve at pseudotime
#'   \code{0}. The weighting curve decreases from \code{1} to \code{0} over the
#'   non-outlying pseudotime values of shared cells (where outliers are defined
#'   by the \code{1.5*IQR} rule). The exact shape of the curve in this region is
#'   controlled by \code{shrink.method}, and can follow the shape of any
#'   standard kernel function's cumulative density curve (or more precisely,
#'   survival curve, since we require a decreasing function). Different choices
#'   of \code{shrink.method} to have no discernable impact on the final curves,
#'   in most cases.
#'
#' @details When \code{reweight = TRUE}, weights for shared cells are based on
#'   the quantiles of their projection distances onto each curve. The
#'   distances are ranked and converted into quantiles between \code{0} and
#'   \code{1}, which are then transformed by \code{1 - q^2}. Each cell's weight
#'   along a given lineage is the ratio of this value to the maximum value for
#'   this cell across all lineages.
#'
#' @references Hastie, T., and Stuetzle, W. (1989). "Principal Curves."
#'   \emph{Journal of the American Statistical Association}, 84:502-516.
#'
#' @references Street, K., et al. (2018). "Slingshot: cell lineage and
#'   pseudotime inference for single-cell transcriptomics." \emph{BMC Genomics},
#'   19:477.
#'
#' @return An object of class \code{\link{PseudotimeOrdering}} containing the
#'   pseudotime estimates and lineage assignment weights in the \code{assays}.
#'   The \code{reducedDim} and \code{clusterLabels} matrices will be stored in
#'   the \code{\link[TrajectoryUtils]{cellData}}. Additionally, the
#'   \code{metadata} slot will contain an \code{\link[igraph]{igraph}} object
#'   named \code{mst}, a list of parameter values named \code{slingParams}, a
#'   list of lineages (ordered sets of clusters) named \code{lineages}, and a
#'   list of \code{\link[princurve]{principal_curve}} objects named
#'   \code{curves}.
#'
#' @examples
#' data("slingshotExample")
#' rd <- slingshotExample$rd
#' cl <- slingshotExample$cl
#' pto <- slingshot(rd, cl, start.clus = '1')
#'
#' # plotting
#' sds <- as.SlingshotDataSet(pto)
#' plot(rd, col = cl, asp = 1)
#' lines(sds, type = 'c', lwd = 3)
#'
#' @export
#'
setMethod(f = "slingshot",
          signature = signature(data = "matrix",
                                clusterLabels = "character"),
          definition = function(data, clusterLabels,
                                reducedDim = NULL,
                                start.clus = NULL, end.clus = NULL,
                                dist.method = "slingshot",
                                use.median = FALSE,
                                omega = FALSE, omega_scale = 1.5,
                                times = NULL,
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                reassign = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                approx_points = NULL,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, distmat=NULL, ...){
              pto <- getLineages(data, clusterLabels, reducedDim = reducedDim,
                                 start.clus = start.clus, end.clus = end.clus,
                                 dist.method = dist.method,
                                 use.median = use.median, omega = omega,
                                 omega_scale = omega_scale, times = times, distmat=distmat)
              pto <- getCurves(pto,
                               shrink = shrink, extend = extend,
                               reweight = reweight, reassign = reassign,
                               thresh = thresh, maxit = maxit,
                               approx_points = approx_points,
                               stretch = stretch, smoother = smoother,
                               shrink.method = shrink.method,
                               allow.breaks = allow.breaks, ...)
              return(pto)
          }
)

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "matrix",
                                clusterLabels = "matrix"),
          definition = function(data, clusterLabels,
                                reducedDim = NULL,
                                start.clus = NULL, end.clus = NULL,
                                dist.method = "slingshot",
                                use.median = FALSE,
                                omega = FALSE, omega_scale = 1.5,
                                times = NULL,
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                reassign = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                approx_points = NULL,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, distmat=NULL, ...){
              pto <- getLineages(data, clusterLabels, reducedDim = reducedDim,
                                 start.clus = start.clus, end.clus = end.clus,
                                 dist.method = dist.method,
                                 use.median = use.median, omega = omega,
                                 omega_scale = omega_scale, times = times, distmat=distmat)
              pto <- getCurves(pto,
                               shrink = shrink, extend = extend,
                               reweight = reweight, reassign = reassign,
                               thresh = thresh, maxit = maxit,
                               approx_points = approx_points,
                               stretch = stretch, smoother = smoother,
                               shrink.method = shrink.method,
                               allow.breaks = allow.breaks, ...)
              return(pto)
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "SlingshotDataSet",
                                clusterLabels = "ANY"),
          definition = function(data, clusterLabels, ...){
            return(slingshot(data = reducedDim(data),
                             clusterLabels = slingClusterLabels(data), ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "data.frame",
                                clusterLabels = "ANY"),
          definition = function(data, clusterLabels, ...){
            RD <- as.matrix(data)
            rownames(RD) <- rownames(data)
            return(slingshot(data = RD,
                             clusterLabels = clusterLabels, ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "matrix",
                                clusterLabels = "numeric"),
          definition = function(data, clusterLabels, ...){
            return(slingshot(data = data,
                             clusterLabels = as.character(clusterLabels), ...))
          })

#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "matrix",
                                clusterLabels = "factor"),
          definition = function(data, clusterLabels, ...){
            return(slingshot(data = data,
                             clusterLabels = as.character(clusterLabels), ...))
          })


#' @rdname slingshot
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "matrix",
                                clusterLabels = "ANY"),
          definition = function(data, clusterLabels, ...){
            if(missing(clusterLabels)){
              message('No cluster labels provided.',
                             ' Continuing with one cluster.')
              clusterLabels <- rep('1', nrow(data))
            }
            if(! any(c(length(clusterLabels), nrow(clusterLabels)) ==
                    nrow(data))){
                stop("clusterLabels must have length or number of rows equal',
                    'to nrow(data).")
            }
            return(slingshot(data = data, clusterLabels = clusterLabels, ...))
          })


#' @rdname slingshot
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SingleCellExperiment reducedDims<-
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "ClusterExperiment"),
          definition = function(data, clusterLabels,
                                reducedDim = NULL,
                                start.clus = NULL, end.clus = NULL,
                                dist.method = "slingshot",
                                use.median = FALSE,
                                omega = FALSE, omega_scale = 1.5,
                                times = NULL,
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                reassign = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                approx_points = NULL,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, distmat=NULL, ...){
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
                  cl <- clusterExperiment::primaryClusterNamed(data)
              }else{
                  if(length(clusterLabels)==1){
                      if(clusterLabels %in% colnames(colData(data))){
                          cl <- colData(data)[[as.character(clusterLabels)]]
                      }else{
                          if(clusterLabels %in% colnames(
                            clusterExperiment::clusterMatrix(data))){
                              cl <- clusterExperiment::clusterMatrixNamed(
                                data)[,as.character(clusterLabels)]
                          }else{
                              stop(clusterLabels,' not found in colData(data)',
                                  ' or clusterMatrix(data).')
                          }
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
              pto <- slingshot(data = rd, clusterLabels = cl,
                               reducedDim = NULL,
                               start.clus = start.clus, end.clus = end.clus,
                               dist.method = dist.method,
                               use.median = use.median, omega = omega,
                               omega_scale = omega_scale, times = times,
                               shrink = shrink, extend = extend,
                               reweight = reweight, reassign = reassign,
                               thresh = thresh, maxit = maxit,
                               stretch = stretch, smoother = smoother,
                               approx_points = approx_points,
                               shrink.method = shrink.method,
                               allow.breaks = allow.breaks, distmat=distmat,  ...)
              # combine slingshot output with SCE
              sce <- data
              colData(sce)$slingshot <- pto
              pst <- slingPseudotime(pto)
              colnames(pst) <- paste0('slingPseudotime_',seq_len(ncol(pst)))
              colData(sce) <- cbind(colData(sce), pst)
              return(sce)
          })



#' @rdname slingshot
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment colData<-
#' @importFrom SingleCellExperiment reducedDims
#' @importFrom SingleCellExperiment reducedDims<-
#' @export
setMethod(f = "slingshot",
          signature = signature(data = "SingleCellExperiment"),
          definition = function(data, clusterLabels,
                                reducedDim = NULL,
                                start.clus = NULL, end.clus = NULL,
                                dist.method = "slingshot",
                                use.median = FALSE,
                                omega = FALSE, omega_scale = 1.5,
                                times = NULL,
                                shrink = TRUE,
                                extend = 'y',
                                reweight = TRUE,
                                reassign = TRUE,
                                thresh = 0.001, maxit = 15, stretch = 2,
                                approx_points = NULL,
                                smoother = 'smooth.spline',
                                shrink.method = 'cosine',
                                allow.breaks = TRUE, distmat=NULL, ...){
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
                  message('No cluster labels provided. Continuing with ',
                          'one cluster.')
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
              pto <- slingshot(data = rd, clusterLabels = cl,
                               reducedDim = NULL,
                               start.clus = start.clus, end.clus = end.clus,
                               dist.method = dist.method,
                               use.median = use.median, omega = omega,
                               omega_scale = omega_scale, times = times,
                               shrink = shrink, extend = extend,
                               reweight = reweight, reassign = reassign,
                               thresh = thresh, maxit = maxit,
                               stretch = stretch, smoother = smoother,
                               approx_points = approx_points,
                               shrink.method = shrink.method,
                               allow.breaks = allow.breaks, distmat=distmat, ...)
              # combine slingshot output with SCE
              sce <- data
              colData(sce)$slingshot <- pto
              pst <- slingPseudotime(pto)
              colnames(pst) <- paste0('slingPseudotime_',seq_len(ncol(pst)))
              colData(sce) <- cbind(colData(sce), pst)
              return(sce)
          })


