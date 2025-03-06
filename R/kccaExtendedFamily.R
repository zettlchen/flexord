#20.01.25

#' Extending K-Centroids Clustering to (Mixed-with-)Ordinal Data
#' 
#' @description
#' This wrapper creates objects of class `"kccaFamily"`,
#' which can be used with [flexclust::kcca()] to conduct K-centroids
#' clustering using the following methods:
#'  - **kModes** (after Weihs et al., 2005)
#'  - **kGower** (Gower's distance after Kaufman & Rousseeuw, 1990,
#'                and a user specified centroid)
#'  - **kGDM2** (GDM2 distance after Walesiak et al., 1993, and a
#'              user specified centroid)
#'              
#' @usage
#' kccaExtendedFamily(which = c('kModes', 'kGDM2', 'kGower'),
#'                    cent = NULL,
#'                    preproc = NULL,
#'                    xrange = NULL,
#'                    xmethods = NULL,
#'                    trim = 0, groupFun = 'minSumClusters')
#'              
#' @details
#' **Wrappers** for defining families are obtained by specifying
#' `which` using:
#'
#' * `which='kModes'` creates an object for **kModes** clustering,
#'   i.e., K-centroids clustering using Simple Matching Distance
#'   (counts of disagreements) and modes as centroids.  Argument
#'   `cent` is ignored for this method.
#' 
#' * `which='kGower'` creates an object for performing clustering
#'   using Gower's method as described in Kaufman & Rousseeuw (1990):
#'
#'   - Numeric and/or ordinal variables are scaled by
#'     \eqn{\frac{\mathbf{x}-\min{\mathbf{x}}}{\max{\mathbf{x}-\min{\mathbf{x}}}}}.
#'     Note that for ordinal variables the internal coding with values
#'     from 1 up to their maximum level is used.
#' 
#'   - Distances are calculated for each column (Euclidean distance,
#'     `distEuclidean`, is recommended for numeric, Manhattan
#'     distance, `distManhattan` for ordinal, Simple Matching
#'     Distance, `distSimMatch` for categorical, and Jaccard distance,
#'     `distJaccard` for asymmetric binary variables), and they are
#'     summed up as:
#'
#'     \deqn{d(x_i, x_k) = \frac{\sum_{j=1}^p \delta_{ikj} d(x_{ij},
#'     x_{kj})}{\sum_{j=1}^p \delta_{ikj}}}
#'
#'     where \eqn{p} is the number of variables and with the weight
#'     \eqn{\delta_{ikj}} being 1 if both values \eqn{x_{ij}} and
#'     \eqn{x_{kj}} are not missing, and in the case of asymmetric
#'     binary variables, at least one of them is not 0.
#'
#'     The columnwise distances used can be influenced in two ways: By
#'     passing a character vector of length \eqn{p} to `xmethods` that
#'     specifies the distance for each column.  Options are:
#'     `distEuclidean`, `distManhattan`, `distJaccard`, and
#'     `distSimMatch`.  Another option is to not specify any methods
#'     within `kccaExtendedFamily`, but rather pass a `"data.frame"`
#'     as argument `x` in `kcca`, where the class of the column is
#'     used to infer the distance measure. `distEuclidean` is used on
#'     numeric and integer columns, `distManhattan` on columns that
#'     are coded as ordered factors, `distSimMatch` is the default for
#'     categorically coded columns, and `distJaccard` is the default
#'     for binary coded columns.
#'   
#'     For this method, if `cent=NULL`, a general purpose optimizer
#'     with `NA` omission is applied for centroid calculation.
#'   
#' * `which='kGDM2'` creates an obejct for clustering using the GDM2
#'   distance for ordinal variables. The GMD2 distance was first
#'   introduced by Walesiak et al. (1993), and adapted in Ernst et
#'   al. (2025), as the distance measure within [flexclust::kcca()].
#' 
#'   This distance respects the ordinal nature of a variable by
#'   conducting only relational operations to compare values, such as
#'   \eqn{\leq}, \eqn{\geq} and \eqn{=}. By obtaining the relative
#'   frequencies and empirical cumulative distributions of \eqn{x}, we
#'   allow for comparison of two arbitrary values, and thus are able
#'   to conduct K-centroids clustering. For more details, see Ernst et
#'   al. (2025).
#' 
#' Also for this method, if `cent=NULL`, a general purpose optimizer
#' with `NA` omission will be applied for centroid calculation.
#' 
#' **Scale handling**.
#' In `'kModes'`, all variables are treated as unordered factors.
#' In `'kGDM2'`, all variables are treated as ordered factors, with strict assumptions
#' regarding their ordinality.
#' `'kGower'` is currently the only method designed to handle mixed-type data. For ordinal
#' variables, the assumptions are more lax than with GDM2 distance.
#' 
#' **NA handling**.
#' NA handling via omission and upweighting non-missing variables is currently
#' only implemented for `'kGower'`. Within `'kModes'`, the omission of NA responses
#' can be avoided by coding missings as valid factor levels. For `'kGDM2'`, currently
#' the only option is to omit missing values completely.
#' 
#' @param which One of either `'kModes'`, `'kGDM2'` or `'kGower'`, the three
#'     predefined methods for K-centroids clustering.  For more
#'     information on each of them, see the Details section.
#' @param cent Function for determining cluster centroids.
#'
#'     This argument is ignored for `which='kModes'`, and `centMode`
#'     is used.  For `'kGDM2'` and `'kGower'`, `cent=NULL` defaults to
#'     a general purpose optimizer.
#' @param preproc Preprocessing function applied to the data before
#'     clustering.
#'
#'     This argument is ignored for `which='kGower'`. In this case,
#'     the default preprocessing proposed by Gower (1971) and Kaufman
#'     & Rousseeuw (1990) is conducted. For `'kGDM2'` and `'kModes'`,
#'     users can specify preprocessing steps here, though this is not
#'     recommended.
#' @param xrange The range of the data in `x`. Options are:
#'
#' - `"all"`: uses the same minimum and maximum value for each column
#'     of `x` by determining the whole range of values in the data
#'     object `x`.
#' 
#' - `"columnwise"`: uses different minimum and maximum values for
#'     each column of `x` by determining the columnwise ranges of
#'     values in the data object `x`.
#'
#' - A vector of `c(min, max)`: specifies the same minimum and maximum
#'     value for each column of `x`.
#' 
#' - A list of vectors `list(c(min1, max1), c(min2, max2),...)` with
#'     length `ncol(x)`: specifies different minimum and maximum
#'     values for each column of `x`. 
#'
#' This argument is ignored for `which='kModes'`. `xrange=NULL`
#' defaults to `"all"` for `'kGDM2'`, and to `"columnwise"` for
#' `'kGower'`.
#' 
#' @param xmethods An optional character vector of length `ncol(x)`
#'     that specifies the distance measure for each column of
#'     `x`. Currently only used for `'kGower'`. For `'kGower'`,
#'     `xmethods=NULL` results in the use of default methods for each
#'     column of `x`. For more information on allowed input values,
#'     and default measures, see the Details section.
#' @param trim Proportion of points trimmed in robust clustering, wee
#'     [flexclust::kccaFamily()].
#' @param groupFun A character string specifying the function for
#'     clustering.
#'
#'     Default is `'minSumClusters'`, see [flexclust::kccaFamily()].
#'   
#' @return An object of class `"kccaFamily"`.
#' 
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).
#'   *Ordinal Clustering with the flex-Scheme.*
#'   Austrian Journal of Statistics. _Submitted manuscript_.
#' - Gower, JC (1971).
#'   *A General Coefficient for Similarity and Some of Its Properties.*
#'   Biometrics, 27(4), 857-871.
#'   \doi{10.2307/2528823}
#' - Kaufman, L, Rousseeuw, P (1990).
#'   *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics.
#'   \doi{10.1002/9780470316801}
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   \doi{10.1016/j.csda.2005.10.006}
#' - Walesiak, M (1993). *Statystyczna Analiza Wielowymiarowa w Badaniach Marketingowych.*
#'   Wydawnictwo Akademii Ekonomicznej, 44-46.
#' - Weihs, C, Ligges, U, Luebke, K, Raabe, N (2005).
#'   *klaR Analyzing German Business Cycles.* In: Data Analysis and
#'   Decision Support, Springer: Berlin. 335-343.
#'   \doi{10.1007/3-540-28397-8_36}
#' 
#' @examples
#' # Example 1: kModes
#' set.seed(123)
#' dat <- data.frame(cont = sample(1:100, 10, replace=TRUE)/10,
#'                   bin_sym = as.logical(sample(0:1, 10, replace=TRUE)),
#'                   bin_asym = as.logical(sample(0:1, 10, replace=TRUE)),                     
#'                   ord_levmis = factor(sample(1:5, 10, replace=TRUE),
#'                                       levels=1:6, ordered=TRUE),
#'                   ord_levfull = factor(sample(1:4, 10, replace=TRUE),
#'                                        levels=1:4, ordered=TRUE),
#'                   nom = factor(sample(letters[1:4], 10, replace=TRUE),
#'                                levels=letters[1:4]))
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kModes'))
#' 
#' # Example 2: kGDM2
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kGDM2',
#'                                                     xrange='columnwise'))
#' # Example 3: kGower
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower'))
#' nas <- sample(c(TRUE,FALSE), prod(dim(dat)), replace=TRUE, prob=c(0.1,0.9)) |> 
#'    matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower',
#'                                            xrange='all'))
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower',
#'                                            xmethods=c('distEuclidean',
#'                                                       'distEuclidean',
#'                                                       'distJaccard',
#'                                                       'distManhattan',
#'                                                       'distManhattan',
#'                                                       'distSimMatch')))
#' #the case where column 2 is a binary variable, but is symmetric
#' 
#' @seealso
#' [flexclust::kcca()],
#' [flexclust::stepFlexclust()],
#' [flexclust::bootFlexclust()]
#'
#' @import flexclust

#' @export
kccaExtendedFamily <- function(which=c('kModes', 'kGDM2', 'kGower'),
                               cent=NULL, preproc=NULL,
                               xrange=NULL, xmethods=NULL,
                               trim=0, groupFun='minSumClusters') { #the last two are unused leftovers from kcca, should probably not provide them
  which <- match.arg(which)
  if(which=='kModes') {
    
    distGen <- NULL
    dstfnc <- distSimMatch
    # was this:
    #cent <- centMode
    # changed to this, so R CMD check doesn't complain.
    cent <- function(x, genDist) {
        centMode(x)
    }

  }
  
  if(which=='kGDM2') {
    
    if(is.null(xrange)) xrange <- 'all'
    
    rng <- .rangeMatrix(xrange)
    
    if(is.null(preproc)) preproc <- function(x, xclass) x #added here because xclass also runs for kGDM2, and thus this preproc also needs its own function environment
    
    distGen <- function(x, xclass) .projectIntofx(x, rangeMatrix=rng)
    dstfnc <- distGDM2
    
    if(is.null(cent)) {
      cent <- function(x, genDist){
        flexclust::centOptim(x, dist = \(y, centers) {
          distGDM2(y, centers, genDist=genDist)
        }) #filler cent, will be recreated in the function
      }
    }
    
  }
  
  if(which=='kGower') {
    
    if(is.null(xrange)) xrange <- 'columnwise'
    
    rng <- .rangeMatrix(xrange)
    
    #old, archived preproc: preproc <- function(x) .ScaleGower(x, rangeMatrix=rng)
    
    if(is.null(xmethods)) {
      # sorry I dont think we should have a warning for the default case
      #warning('No column-wise distance measures specified, default measures
      #      for each column will be used.')
      distGen <- function(x, xclass) {
        #I apologize for the use of parent.frame(), but didn't know how else to fix it.
        #however I do think it's ok here because 1) 'xclass' is not a generic method
        #and 2) because I only use them once in the beginning (unlike the dists),
        #so I can't really get lost in the frames
        #xcls <- get('xclass', parent.frame())
        .ChooseVarDists(xclass)
      }
      preproc <- function(x, xclass) {
        #xcls <- get('xclass', parent.frame())
        .ScaleVarSpecific(x, rangeMatrix=rng,
                                    xclass=xclass)
      }
    } else {
      distGen <- function(x, xclass) {
        if(!all(xmethods %in% c('distEuclidean', 'distManhattan',
                                'distSimMatch', 'distJaccard')))
          stop('Specified columnwise xmethod not implemented!')
        return(xmethods)
      }
      preproc <- function(x, xclass) .ScaleVarSpecific(x, rangeMatrix=rng,
                                               xclass=xmethods)
    }
    
    dstfnc <- distGower
    
    #the default combo in the paper for distGower was centMin. now that x is scaled in the beginning, I think centOptim is sufficient
    if(is.null(cent)) {
      cent <- function(x, genDist){
        centOptimNA(x, dist = \(y, centers) {
          distGower(y, centers, genDist=genDist)
        }) #filler cent, will be recreated in the function
      }
    }
    
  }
 

  newgendist <- function(x, family) {
    kccaExtendedFamilyGenDist(x, family, genDist = distGen)
  }



  flexclust::kccaFamily(name=which,
                        dist=dstfnc,
                        cent=cent,
                        genDist=newgendist,
                        preproc=preproc,
                        trim=trim, groupFun=groupFun)
                              
}



#' Recreate the Family Object with Updated Distance, Centroid, ... Functions 
#' @param x the data set that kcca was called with
#' @param family the original family
#' @param genDist a function that generates a vector of distance functions
#' @noRd
kccaExtendedFamilyGenDist = function(x, family, genDist) {
  if(is.data.frame(x)) {
    xclass <- sapply(x, data.class)
  } else {
    xclass <- rep('numeric', ncol(x))
  }

  origDist <- family@dist
  origCent <- family@cent
  origPreproc <- family@preproc

  newpreproc <- if("xclass" %in% names(formals(origPreproc))) {
    function(x) origPreproc(x, xclass = xclass)
  } else {
    origPreproc
  }

  newdist <- if("genDist" %in% names(formals(origDist))) {
    function(x, centers) origDist(x, centers, genDist = generated_dist)
  } else {
    origDist
  }

  newcent <- if("genDist" %in% names(formals(origCent))) {
    function(x) origCent(x, genDist = generated_dist)
  } else {
    origCent
  }

  family_new <- kccaFamily(
    name     = family@name,
    dist     = newdist,
    cent     = newcent,
    genDist = family@genDist,
    preproc  = newpreproc,
    trim     = family@trim,
    groupFun = family@groupFun)

  family <- family_new

  x <- data.matrix(x) #previously: x <- as(x, "matrix")  
  x <- family@preproc(x)

  if (!is.null(genDist)) {
      generated_dist <-
          if("xclass" %in% names(formals(genDist))) {
              genDist(x, xclass = xclass)
          } else {
              genDist(x)
          }
  }

  family
}


