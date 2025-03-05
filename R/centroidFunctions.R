#' @rdname centroids
#' @export
centMode <- function(x) {
  apply(x, 2, \(cat) {
    matrix(tabulate(cat), nrow = 1) |>
      max.col(ties.method = 'random') # Other options: first, last
  })
}

#' @rdname centroids
#' @export
centMin <- function(x, dist, xrange=NULL) {
  if(is.null(xrange)) xrange <- 'all'
  
  rng <- .rangeMatrix(xrange=xrange)(x)
  
  lvls <- apply(rng, 2, \(y) matrix(seq(y[1], y[2])),
                simplify = F)
  cntrs <- vector(length = ncol(x)) |> 
    stats::setNames(colnames(x))
  for (i in 1:ncol(x)) {
    d <- dist(x[, i, drop=F], lvls[[i]]) |> 
      colSums() |> matrix(nrow=1)
    cntrs[i] <- max.col(-d, ties.method='random')
  }
  return(cntrs)
}

#' @rdname centroids
#' @export
centOptimNA <- function(x, dist) {
  foo <- function(p)
    sum(dist(x, matrix(p, nrow=1)), na.rm=TRUE)
  
  stats::optim(colMeans(x, na.rm=TRUE), foo)$par
}

#' Centroid Functions for K-Centroids Clustering of (Ordinal) Categorical/Mixed Data
#'
#' @description
#' Functions to calculate cluster centroids for K-centroids clustering that extend the
#' options available in package \pkg{flexclust}.
#' 
#' `centMode` calculates centroids based on the mode of each variable.
#' `centMin` minimizes the applied distance metric to find centroids within a specified range.
#' `centOptimNA` replicates the exact behaviour of [flexclust::centOptim()], just with NA removal.
#' 
#' These functions are designed for use within [flexclust::kcca()] or functions that are built
#' upon it. Their use is easiest via the wrapper `kccaExtendedFamily`, for more information
#' see [kccaExtendedFamily()].
#'
#' @details
#' - **`centMode`**: Column-wise modes are used as centroids, and ties are
#'  broken randomly. In combination with Simple Matching Distance (`distSimMatch`),
#'  this results in the `kmodes` algorithm.
#'
#' - **`centMin`**: Column-wise centroids are calculated by minimizing the
#'  specified distance measure between a value of `x`, and all possible levels of `x`.
#'  
#'  - **centOptimNA**: Column-wise centroids are calculated by minimizing the
#'  specified distance measure via a general purpose optimizer. Unlike in [flexclust::centOptim()],
#'  NAs are removed from the starting search values.
#'
#' @param x A numeric matrix or data frame. Categorical/ordinal variables
#'    need to be coded as `1:length(levels(x[,i]))` in steps of one.
#' @param dist The distance measure function to be minimized in `centMin`.
#' @param xrange Range specification for the variables. Currently only used for `centMin`.
#'   Possible values are:
#'   - `NULL` (default): defaults to `'all'`.
#'   - `'all'`: Uses the range of the entire data set.
#'   - `'columnwise'`: Uses column-specific ranges.
#'   - A numeric vector of c(min, max): Applies the specified range to all columns.
#'   - A list of numeric vectors of c(min,max): Uses the user-specified columnwise
#'       ranges, the length of the list must be equal to the number of columns to be scaled.
#'
#' @return 
#' A named numeric vector containing the centroid values for each column of `x`.
#' 
#' @seealso
#' [kccaExtendedFamily()],
#' [flexclust::kcca()]
#'
#' @importFrom stats setNames optim
#'
#' @examples
#' # Example: Mode as centroid
#' dat <- data.frame(A = rep(2:5, 2),
#'                   B = rep(1:4, 2),
#'                   C = rep(c(1, 2, 4, 5), 2))
#' centMode(dat)
#' ## within kcca
#' flexclust::kcca(dat, 3, family=kccaExtendedFamily('kModes')) #default centroid
#' 
#' # Example: Centroid is level for which distance is minimal
#' centMin(dat, flexclust::distManhattan, xrange = 'all')
#' ## within kcca
#' flexclust::kcca(dat, 3,
#'                 family=flexclust::kccaFamily(dist=flexclust::distManhattan,
#'                                              cent=\(y) centMin(y, flexclust::distManhattan,
#'                                                                xrange='all')))
#'                              
#' # Example: Centroid calculated by general purpose optimizer with NA removal
#' nas <- sample(c(TRUE, FALSE), prod(dim(dat)),
#'               replace=TRUE, prob=c(0.1,0.9)) |> 
#'        matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' centOptimNA(dat, flexclust::distManhattan)
#' ## within kcca
#' flexclust::kcca(dat, 3, family=kccaExtendedFamily('kGower')) #default centroid
#' @name centroids
NULL

