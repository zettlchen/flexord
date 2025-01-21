#21.01.25
#This script contains the code for distSimMatch, and the main documentation
#block for distance_functions.Rd, which is the documentation page for
#all three help pages together.
#distGDM2.R contains code for distGDM2 and .projectIntofx
#distGower.R contains code for distGower, .ChooseVarDists, .ScaleVarSpecific,
# .delta, .distGower_mixedType and .distGower_SingleTypeNoNAs.

#' New Distance Functions for K-centroids clustering of (ordinal) categorical/mixed data
#'
#' @description
#' Functions to calculate the distance between a matrix `x` and an arbitrary
#' matrix `c`, which can be used for K-centroids clustering via `flexclust::kcca`.
#' 
#' `distSimMatch` implements Simple Matching Distance (most frequently
#' used for categorical, or symmetric binary data) into K-centroids
#' clustering.
#' 
#' `distGower` implements Gower's Distance after Gower (1971) and
#' Kaufman\&Rousseeuw (1990) for mixed-type data with missings into K-centroids
#' clustering.
#' 
#' `distGDM2` implements GDM2 distance for ordinal data introduced by
#' Walesiak et\~al. (1993) and adapted to K-centroids clustering by
#' Ernst et\~al. (2025).
#' 
#' These functions are designed for use within `flexclust::kcca` or functions that are built
#' upon it. Their use is easiest via the wrapper `kccaExtendedFamily`, for more information
#' see `?kccaExtendedFamily`. However, they can easily be extended to result
#' in a distance matrix of `x`, see Examples.
#' 
#' @details
#' - `distSimMatch`: Simple Matching Distance between two individuals is calculated as
#'                  the proportion of disagreements out of all variables. Described f.i.
#'                  in Kaufman & Rousseeuw (1990), p.24.
#'                  If this is used in K-centroids analysis in combination with the
#'                  mode (as implemented in `centMode`), this results in the `kModes`
#'                  algorithm.
#'                  A wrapper for this algorithm is found in `kccaExtendedFamily(which='kModes')`.
#' - `distGower`: Distances are calculated for each column \eqn{p} (squared Euclidean distance, `distEuclidean`,
#'   is recommended for numeric, Manhattan distance, `distManhattan` for ordinal, Simple Matching
#'   Distance, `distSimMatch` for categorical, and Jaccard distance, `distJaccard` for asymmetric
#'   binary variables), and they are summed up as:
#'   \deqn{d(x_i, x_k) = \frac{\sum_{j=1}^p \delta_{ikj} d(x_{ij}, x_{kj})}{\sum_{j=1}^p \delta_{ikj}}}
#'   with the weight \eqn{\delta_{ikj}} being 1 if both values \eqn{x_{ij}} and \eqn{x_{kj}} are
#'   not missing, and in the case of asymmetric binary variables, at least one of them is not 0.
#'   Please note that for calculating Gower's distance, scaling of numeric/ordered
#'   variables is required (as f.i. by `.ScaleVarSpecific`).
#'   A wrapper for K-centroids analysis using Gower's distance in combination with a numerically
#'   optimized centroid is found in `kccaExtendedFamily(which='kGower')`.
#' - `distGDM2`: GDM2 distance for ordinal variables conducts only relational operations
#'    on the variables, such as \eqn{\leq}, \eqn{\geq} and \eqn{=}. By translating \eqn{x}
#'    to its relative frequencies and empirical cumulative distributions, we are able to
#'    extend this principle to compare two arbitrary values, and thus use it within K-Centroids
#'    clustering. For more details, see Ernst et\~al. (2025).
#'    A wrapper for this algorithm in combination with a numerically optimized centroid
#'    is found in `kccaExtendedFamily(which='kGDM2')`.
#'    
#' The distances functions presented here can also be used in clustering algorithms that
#' rely on distance matrices (such as hierarchical clustering and PAM), if applied
#' accordingly, see Examples.
#'   
#'
#' @param x A numeric matrix or data frame. Categorical/ordinal variables
#'    need to be coded as `1:length(levels(x[,i]))` in steps of one.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#' `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' @param genDist Additional information on x required for distance calculation.
#' Filled automatically if used within `flexclust::kcca`.
#'    - For `distGower`: A character vector of variable specific distances to be used,
#'                  f.i. as derived from `.ChooseVarDists(x)`. Length
#'                  needs to be equal to `ncol(x)`. Can contain the options:
#'                  - `distEuclidean`: squared Euclidean distance between the
#'                                     scaled variables
#'                  - `distManhattan`: absolute distance between the scaled variables
#'                  - `distJaccard`: counts of zero if both binary variables are
#'                                   equal to 1, and 1 otherwise
#'                  - `distSimMatch`: Simple Matching Distance, i.e. the number of agreements
#'                                    between variables.
#'    - For `distGDM2`: Function for creating a distance function that will be
#'        primed on `x`, such as f.i. `.projectIntofx`. `.projectIntofx` is a helper
#'        function factory that creates a function that will project an arbitrary object
#'        `centers` into the space of the empirical PDF and CDF of `x`.
#'    - For `distSimMatch`: not used.
#' @param xrange Range specification for the variables. Currently only used for `distGDM2`
#'  (as `distGower` expects `x` to be already scaled). Possible values are:
#'   - `NULL` (default): defaults to `'all'`.
#'   - `'all'`: Uses the range of the entire data set.
#'   - `'columnwise'`: Uses column-specific ranges.
#'   - A numeric vector of c(min, max): Applies the specified range to all columns.
#'   - A list of numeric vectors of c(min,max): Uses the user-specified columnwise
#'       ranges, the length of the list must be equal to the number of columns to be scaled.
#'
#' @return
#' A matrix of dimensions `c(nrow(x),nrow(centers))` that contains the distance
#' between each row of `x` from each row of `centers`.
#'
#' @examples
#' # Example 1: Simple Matching Distance
#' set.seed(123)
#' dat <- data.frame(question1 = factor(sample(LETTERS[1:4], 10, replace=TRUE)),
#'                   question2 = factor(sample(letters[1:6], 10, replace=TRUE)),
#'                   question3 = factor(sample(LETTERS[1:4], 10, replace=TRUE)),
#'                   question4 = factor(sample(LETTERS[1:5], 10, replace=TRUE)),
#'                   state = factor(sample(state.name[1:10], 10, replace=TRUE)),
#'                   gender = factor(sample(c('M', 'F', 'N'), 10, replace=TRUE,
#'                                          prob=c(0.45, 0.45, 0.1))))
#' datmat <- data.matrix(dat)
#' initcenters <- datmat[sample(1:10, 3),]
#' distSimMatch(datmat, initcenters)
#' ## within kcca
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kModes'))
#' ## as a distance matrix
#' as.dist(distSimMatch(datmat, datmat))
#' 
#' # Example 2: GDM2 distance
#' distGDM2(datmat, initcenters, genDist=flexord:::.projectIntofx)
#' ## within kcca
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kGDM2'))
#' ## as a distance matrix
#' as.dist(distGDM2(datmat, datmat, genDist=flexord:::.projectIntofx))
#' 
#' # Example 3: Gower's distance
#' # Ex. 3.1: single variable type case with no missings:
#' xcls <- flexord:::.ChooseVarDists(datmat) #all Euclidean (on dat, it would default to all Simple Matching)
#' datscld <- flexord:::.ScaleVarSpecific(datmat, xclass=xcls,
#'                                        xrange=list(c(1,4), c(1,6), c(1,4),
#'                                                    c(1,5), c(1,10), c(1,3)))
#' initcentscld <- datscld[sample(1:10, 3),]
#' distGower(datscld, initcentscld, genDist=xcls)
#' ## within kcca
#' flexclust::kcca(datmat, 3, kccaExtendedFamily('kGower')) #turns into kmeans with scaling
#' 
#' # Ex. 3.2: single variable type case with missing values:
#' nas <- sample(c(TRUE,FALSE), prod(dim(dat)), replace=TRUE, prob=c(0.1,0.9)) |> 
#'    matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' #repeat the steps from above...or just do:
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower', cent=centMode)) #turns into kModes with upweighting of present values
#' 
#' #Ex. 3.3: mixed variable types (with or without missings): 
#' dat <- data.frame(cont = sample(1:100, 10, replace=TRUE)/10,
#'                   bin_sym = as.logical(sample(0:1, 10, replace=TRUE)),
#'                   bin_asym = as.logical(sample(0:1, 10, replace=TRUE)),                     
#'                   ord_levmis = factor(sample(1:5, 10, replace=TRUE),
#'                                       levels=1:6, ordered=TRUE),
#'                   ord_levfull = factor(sample(1:4, 10, replace=TRUE),
#'                                        levels=1:4, ordered=TRUE),
#'                   nom = factor(sample(letters[1:4], 10, replace=TRUE),
#'                                levels=letters[1:4]))
#' dat[nas] <- NA
#' xcls <- flexord:::.ChooseVarDists(dat)
#' datmat <- flexord:::.ScaleVarSpecific(data.matrix(dat), xclass=xcls,
#'                                       xrange='columnwise')
#' initcenters <- datmat[sample(1:10, 3),]
#' distGower(datmat, initcenters, genDist=xcls)                  
#' ## within kcca
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower'))
#' ## as a distance matrix
#' distGower(datmat, datmat, genDist=xcls) |> as.dist()
#' ## as a distance matrix
#' 
#' @seealso
#' - [flexclust::kcca()](https://cran.r-project.org/package=flexclust)
#' - [klaR::kmodes()](https://cran.r-project.org/package=klaR)
#' - [cluster::daisy()](https://search.r-project.org/CRAN/refmans/cluster/html/daisy.html)
#' - [clusterSim::dist.GDM()](https://search.r-project.org/CRAN/refmans/clusterSim/html/dist.GDM.html)
#'
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Gower, JC (1971).
#'   *A General Coefficient for Similarity and Some of Its Properties.*
#'   Biometrics, 27(4), 857-871.
#'   [doi:10.2307/2528823](https://doi.org/10.2307/2528823)
#' - Kaufman, L, Rousseeuw, P (1990).
#'   *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics.
#'   [doi:10.1002/9780470316801](https://doi.org/10.1002/9780470316801)
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   [doi:10.1016/j.csda.2005.10.006](https://doi.org/10.1016/j.csda.2005.10.006)
#' - Kaufman L, Rousseeuw, P (1990.) *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics, New York: John Wiley \& Sons.
#'   [doi:10.1002/9780470316801](https://doi.org/10.1002/9780470316801)
#' - Walesiak, M (1993). *Statystyczna Analiza Wielowymiarowa w Badaniach Marketingowych.*
#'   Wydawnictwo Akademii Ekonomicznej, 44-46.
#' - Weihs, C, Ligges, U, Luebke, K, Raabe, N (2005). *klaR Analyzing German Business Cycles*.
#'   In Baier D, Decker, R, Schmidt-Thieme, L (eds.). Data Analysis and Decision Support,
#'   335-343.Berlin: Springer-Verlag.
#'   [doi:10.1007/3-540-28397-8_36](https://doi.org/10.1007/3-540-28397-8_36)
#' @name distance_functions
NULL

#' @rdname distance_functions
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

