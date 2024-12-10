#06.12.24
#Other, non-refined centroid methods were written: `centOptimLB`, `centMedoid`.  
# TODO: adapt examples when distGower and the wrapper are done. should prolly be
#       in the wrapper anyways

#' Additional Centroid Functions for K-centroids clustering of (ordinal) categorical data
#'
#' @usage
#' centMode(x)
#' centMin(x, dist, xrange=NULL)
#'
#' @description
#' Functions to calculate cluster centroids that extend the options available
#' in `flexclust`.
#' 
#' `centMode` calculates centroids based on the mode of each variable.
#' `centMin` minimizes the applied distance metric to find centroids within a specified range.
#' 
#' These functions are designed for use within `flexclust::kcca` or functions that are built
#' upon it.
#'
#' @details
#' - **`centMode`**: Column-wise modes are used as centroids, and ties are
#'  broken randomly. In combination with Simple Matching Distance (`distSimMatch`),
#'  this results in the `kmodes` algorithm.
#'
#' - **`centMin`**: Column-wise centroids are calculated by minimizing the
#'  specified distance measure between a value of `x`, and all possible levels of `x`.
#'
#' @param x A numeric matrix or data frame. Categorical/ordinal variables
#'    need to be coded as `1:length(levels(x[,i]))`.
#' @param dist The distance measure function to be minimized in `centMin`.
#' @param xrange Range specification for the variables. Possible values are:
#'   - `NULL` (default): defaults to `'data range'`.
#'   - `'data range'`: Uses the range of the entire dataset.
#'   - `'variable specific'`: Uses column-specific ranges.
#'   - A numeric vector: Applies the same range to all columns.
#'   - A list of numeric vectors: Specifies column-specific ranges.
#'
#' @return 
#' A named numeric vector containing the centroid values for each variable.
#' 
#' @seealso - [flexclust::kcca()](https://cran.r-project.org/package=flexclust)
#'
#' @examples
#' # Example: Mode as centroid
#' dat <- data.frame(A = rep(2:5, 2),
#'                   B = rep(1:4, 2),
#'                   C = rep(c(1, 2, 4, 5), 2))
#' centMode(dat)
#' ## within kcca
#' kcca(x, 3, family=kccaFamily(dist=distSimMatch,
#'                              cent=centMode))
#' 
#' # Example: Centroid is level for which distance is minimal
#' centMin(y, distGower, xrange = 'data range')
#' ## within kcca
#'kcca(x, 3, family=kccaFamily(dist= distGower,
#'                              cent=\(y) centMin(y, distGower))
#'
#' @export
centMode <- function(x) {
  apply(x, 2, \(cat) {
    matrix(tabulate(cat), nrow = 1) |>
      max.col(ties.method = 'random') # Other options: first, last
  })
}

#helper function for the different options regarding the xrange parameter
.rangeMatrix <- function(xrange) {
  if(all(xrange=='data range')) {
    rng <- function(x) {
      rep(range(x), ncol(x)) |>
        matrix(nrow=2)
    }
  } else if(all(xrange=='variable specific')) {
    rng <- function(x) {
      apply(x, 2, range)
    }
  } else if(is.vector(xrange, mode='numeric')) {
    rng <- function(x) {
      rep(xrange, ncol(x)) |>
        matrix(nrow=2)
    }
  } else {
    rng <- function(x) {
      if(length(xrange) != ncol(x))
        stop('Either supply 1 range vector, or list of ranges for all variables')
      unlist(xrange) |> matrix(nrow=2)
    }
  }
  return(rng)
}

#' @export
centMin <- function(x, dist, xrange=NULL) {
  if(is.null(xrange)) xrange <- 'data range'
  
  rng <- .rangeMatrix(xrange=xrange)(x)
  
  lvls <- apply(rng, 2, \(y) matrix(seq(y[1], y[2])),
                simplify = F)
  cntrs <- vector(length = ncol(x)) |> 
    setNames(colnames(x))
  for (i in 1:ncol(x)) {
    cntrs[i] <- dist(x[, i, drop = FALSE], lvls[[i]][, 1, drop = FALSE],
                     xrange = rng[, i]) |>
      colSums() |> setNames(lvls[[i]]) |> 
      which.min() |> #automatically picks first
      names() |> as.numeric()
  }
  return(cntrs)
}
