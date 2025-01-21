#09.12.24
#Simple Matching Distance
#have decided to keep the distance and the centroid in separate scripts
#maybe together in docu?
#Other option: the docu for the centroids remains together

#' Simple Matching Distance and Kmodes for K-centroids clustering
#'
#'@usage
#' distSimMatch(x, centers)
#' kccaFamilyKmodes(preproc=NULL, trim=0, groupFun='minSumClusters')
#'
#'@description
#' `distSimMatch` implements Simple Matching Distance (most frequently
#' used for categorical, or symmetric binary data) into K-centroids
#' clustering via package `flexclust`.
#' 
#' Together with modes as centroids (implemented in `centMode`), this
#' results in the `kmodes`-algorithm. Implementations
#' are readily available (see f.i. package `klaR`), but adapting it into
#' the flex-scheme facilitates benchmarking between methods.
#' `kccaFamilyKmodes` is a wrapper for the `kmodes` clustering family
#' to be used within `flexclust::kcca`.
#' 
#' @details
#' Simple Matching Distance between two individuals is calculated as
#' the proportion of disagreements out of all variables. Described f.i.
#' in Kaufman & Rousseeuw (1990), p.24.
#' `centMode` calculates the column-wise most frequent values of a dataset.
#'
#' @param x A matrix of categorical (or binary) data points. They can be
#' coded as character or as numeric. If the latter is chosen, they need
#' to be coded as `1:length(levels(x[,i]))`. Note that for use within
#' `flexclust::kcca`, numeric coding is recommended, as a character matrix
#' will be converted to numeric via `data.matrix()`.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#' `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' Within `kcca`, this is filled with the cluster centers.
#' @param preproc See [flexclust::kccaFamily()]: Preprocessing function to be
#' applied to the data before clustering.
#' @param trim See [flexclust::kccaFamily()]: Proportion of points trimmed in
#' robust clustering.
#' @param groupFun See [flexclust::kccaFamily()]: A character string specifying
#' the function for clustering.
#'   Default is `'minSumClusters'`.
#'
#' @return
#' - `distSimMatch`: A distance matrix for each row in `x` from each row in `centers`
#'    with dimensions `c(nrow(x), nrow(centers)`.
#' - `kccaFamilyKmodes`: A custom `kccaFamily` for `kmodes` clustering, i.e.
#'    K-centroids clustering using Simple Matching Distance and modes for centroids.
#'    To be used within `flexclust::kcca` or functions that build upon it.
#'
#' @examples
#' ## Example usage
#' #creating a distance matrix between two matrices based on GDM2 distance
#' dat <- data.frame(option1 = sample(LETTERS[1:4], 10, replace=TRUE),
#'                   option2 = sample(letters[1:6], 10, replace=TRUE),
#'                   state = sample(state.name[1:10], 10, replace=TRUE))
#' datmat <- data.matrix(dat)
#' initcenters <- datmat[sample(1:10, 3),]
#' distSimMatch(datmat, initcenters)
#' 
#' #a classical distance matrix for a single matrix can be obtained via:
#' as.dist(distSimMatch(datmat, datmat))
#' 
#' #K-centroids clustering via kmodes
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kModes'))
#' 
#' @seealso
#' - [klaR::kmodes()](https://cran.r-project.org/package=klaR)
#' - [flexclust::kcca()](https://cran.r-project.org/package=flexclust)
#'
#' @references
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   [doi:10.1016/j.csda.2005.10.006](https://doi.org/10.1016/j.csda.2005.10.006)
#' - Kaufman L, Rousseeuw, P (1990.) *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics, New York: John Wiley \& Sons.
#'   [doi:10.1002/9780470316801](https://doi.org/10.1002/9780470316801)
#' - Weihs, C, Ligges, U, Luebke, K, Raabe, N (2005). *klaR Analyzing German Business Cycles*.
#'   In Baier D, Decker, R, Schmidt-Thieme, L (eds.). Data Analysis and Decision Support,
#'   335-343.Berlin: Springer-Verlag.
#'   [doi:10.1007/3-540-28397-8_36](https://doi.org/10.1007/3-540-28397-8_36)
#'
#' @export
distSimMatch <- function (x, centers) {
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for(k in 1:nrow(centers)) {
    z[,k] <- colMeans(t(x) != centers[k,])
  }
  z
}

#still creating a separate wrapper here as well, then it's easier for me
#afterward to all smush them into one
#' @export
kccaFamilyKmodes <- function(preproc=NULL, trim=0, groupFun='minSumClusters') {
  flexclust::kccaFamily(name='kmodes',
                        dist=distSimMatch,
                        cent=centMode,
                        genDist=NULL,
                        preproc=preproc, trim=trim, groupFun=groupFun)
}
