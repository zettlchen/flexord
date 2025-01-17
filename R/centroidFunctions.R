#Stand: 17.01.25

#Other centroid methods were drafted: `centOptimLB`, `centMedoid`.  
# TODO: adapt examples when distGower and the wrapper are done. should prolly be
#       in the wrapper anyways

#' Additional Centroid Functions for K-centroids clustering of (ordinal) categorical data
#'
#' @usage
#' centMode(x)
#' centMin(x, dist, xrange=NULL)
#' centOptimNA(x, dist)
#'
#' @description
#' Functions to calculate cluster centroids that extend the options available
#' in `flexclust`.
#' 
#' `centMode` calculates centroids based on the mode of each variable.
#' `centMin` minimizes the applied distance metric to find centroids within a specified range.
#' `centOptimNA` replicates the exact behaviour of `flexclust::centOptim`, just with NA removal.
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
#'  - **centOptimNA**: Column-wise centroids are calculated by minimizing the
#'  specified distance measure via a general purpose optimizer. Unlike in `flexclust::centOptim`,
#'  NAs are removed from the starting search values.
#'
#' @param x A numeric matrix or data frame. Categorical/ordinal variables
#'    need to be coded as `1:length(levels(x[,i]))` in steps of one.
#' @param dist The distance measure function to be minimized in `centMin`.
#' @param xrange Range specification for the variables. Possible values are:
#'   - `NULL` (default): defaults to `'all'`.
#'   - `'all'`: Uses the range of the entire data set.
#'   - `'columnwise'`: Uses column-specific ranges.
#'   - A numeric vector of c(min, max): Applies the specified range to all columns.
#'   - A list of numeric vectors of c(min,max): Uses the user-specified columnwise
#'       ranges, the length of the list must be equal to the number of columns to be scaled.
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
#' kcca(x, 3, family=kccaFamily(dist= distGower,
#'                              cent=\(y) centMin(y, distGower))
#'                              
#' # Example: General purpose optimizer with NA removal
#' nas <- sample(c(T, F), prod(dim(dat)),
#'               replace=T, prob=c(0.1,0.9)) |> 
#'        matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' centOptimNA(dat, distManhattan)
#' @export
centMode <- function(x) {
  apply(x, 2, \(cat) {
    matrix(tabulate(cat), nrow = 1) |>
      max.col(ties.method = 'random') # Other options: first, last
  })
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

#' @export
centOptimNA <- function(x, dist) {
  foo <- function(p)
    sum(dist(x, matrix(p, nrow=1)), na.rm=TRUE)
  
  optim(colMeans(x, na.rm=TRUE), foo)$par
}
