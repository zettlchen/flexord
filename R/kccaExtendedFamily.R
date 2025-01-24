#20.01.25

#' Extending K-Centroids clustering to (mixed-with-)ordinal data
#' 
#' @description
#' This wrapper creates objects of class `flexclust::kccaFamily`,
#' which can be used within `flexclust::kcca` to conduct K-centroids
#' clustering using the following methods:
#'  - **kModes** (after Weihs et\~al., 2005)
#'  - **kGower** (Gower's distance after Kaufman \& Rousseeuw, 1990,
#'                and a user specified centroid)
#'  - **kGDM2** (GDM2 distance after Walesiak et\~al., 1993, and a
#'              user specified centroid)
#'              
#' @usage
#' kccaExtendedFamily(which=c('kModes', 'kGDM2', 'kGower'),
#'                    cent=NULL,
#'                    preproc=NULL,
#'                    xrange=NULL,
#'                    xmethods=NULL,
#'                    trim=0, groupFun='minSumClusters')
#' flexclust::kcca(x, k,
#'                 family=kccaExtendedFamily(which=c('kModes', 'kGDM2', 'kGower'), ...))
#'              
#' @details
#' `kccaExtendedFamily(which='kModes')` creates an object that will perform
#' **kModes** clustering, that is to say K-centroids clustering using Simple
#' Matching Distance (counts of disagreements) and modes as centroids.
#' User-specified alternative centroids in parameter `cent` are ignored for
#' this method.
#' 
#' An object created with `which='kGower` will use Gower's method as described
#' in Kaufman\&Rousseeuw (1990) for clustering:
#' - Numeric and/or ordinal variables are scaled by \eqn{\frac{\mathbf{x}-\min{\mathbf{x}}}{\max{\mathbf{x}-\min{\mathbf{x}}}}}
#'   (please note that ordinal variables have to be coded so that they represent their ranks
#'   in the way that they start with 1 and go in steps of one up to their maximum level \eqn{m}.)
#' - Distances are calculated for each column \eqn{p} (squared Euclidean distance, `distEuclidean`,
#'   is recommended for numeric, Manhattan distance, `distManhattan` for ordinal, Simple Matching
#'   Distance, `distSimMatch` for categorical, and Jaccard distance, `distJaccard` for asymmetric
#'   binary variables), and they are summed up as:
#'   \deqn{d(x_i, x_k) = \frac{\sum_{j=1}^p \delta_{ikj} d(x_{ij}, x_{kj})}{\sum_{j=1}^p \delta_{ikj}}}
#'   with the weight \eqn{\delta_{ikj}} being 1 if both values \eqn{x_{ij}} and \eqn{x_{kj}} are
#'   not missing, and in the case of asymmetric binary variables, at least one of them is not 0.
#'
#'   The columnwise distances used can be influenced twofold: By passing a character
#'   vector of length \eqn{p} to `xmethods` that specifies the distance for each column.
#'   Options are: `distEuclidean`, `distManhattan`, `distJaccard`, and `distSimMatch`.
#'   Another option is to not specify any methods within `kccaExtendedFamily`, but rather
#'   in `kcca` pass a dataframe to `x`, where each column is coded to achieve the desired
#'   treatment. `distEuclidean` is used on numeric and integer columns, `distManhattan`
#'   on columns that are coded as ordered factors, `distSimMatch` is the default for categorically
#'   coded columns, and `distJaccard` is the default for binary coded columns.
#'   
#'   For this method, if `cent=NULL`, a general purpose optimizer with `NA` omission
#'   will be applied for centroid calculation.
#'   
#' An object created with `which=kGDM2` will use the GDM2 distance for ordinal variables,
#' which was first introduced by Walesiak et\~al. (1993), and adapted in Ernst et\~al. (2025),
#' as the distance measure within `flexclust::kcca`.
#' 
#' The principle behind it is that the ordinality of a variable will be
#' respected by conducting only relational operations on them, such as \eqn{\leq}, \eqn{\geq} and \eqn{=}.
#' By translating \eqn{x} to its relative frequencies and empirical cumulative
#' distributions, we are able to extend this principle to compare two arbitrary
#' values, and thus conduct K-Centroids clustering. For more details, see Ernst et\~al. (2025).
#' 
#' Also for this method, if `cent=NULL`, a general purpose optimizer with `NA` omission
#' will be applied for centroid calculation.
#' 
#' **Scale handling**
#' In `kModes`, all variables are treated as unordered factors.
#' In `kGDM2`, all variables are treated as ordered factors, with strict assumptions
#' regarding their ordinality.
#' `kGower` is currently the only method designed to handle mixed-type data. For ordinal
#' variables, the assumptions are more lax than with GDM2 distance.
#' 
#' **NA handling**
#' NA handling via omission and upweighting non-missing variables is currently
#' only implemented for `kGower`. Within `kModes`, the omission of NA responses
#' can be avoided by coding missings as valid factor levels. For `kGDM2`, currently
#' the only option is to omit missing values completely.
#' 
#' @param which One of either `kModes`, `kGDM2` or `kGower`, the three predefined
#'              methods for K-centroids clustering.
#'              For more information on each of them, see the Details section.
#' @param cent Function for determining cluster
#'             centroids (See: [flexclust::kccaFamily()](https://search.r-project.org/CRAN/refmans/flexclust/html/kcca.html)).
#'             This parameter is ignored for `which='kModes'`, and `centMode` is used.
#'             For `kGDM2` and `kGower`, `cent=NULL` defaults to a general purpose optimizer.
#' @param preproc Preprocessing function to be applied to the data before clustering
#'                (See: [flexclust::kccaFamily()](https://search.r-project.org/CRAN/refmans/flexclust/html/kcca.html)).
#'                This parameter is ignored for `which='kGower'`, instead, the default
#'                preprocessing proposed by Gower (1971) and Kaufman\&Rousseeuw (1990)
#'                is conducted. For `kGDM2` and `kModes`, users can specify preprocessing
#'                steps here, though this is not recommended.
#' @param xrange The range of the data in `x`. Options are:
#'               - 'all', which will use the whole range of the data object `x`,
#'               - 'columnwise', which will use columnwise ranges of the data object `x`,
#'               - a vector of `c(min,max)`, which allows the user to specify the
#'                 absolute possible minimum and maximum of the data object `x`,
#'               - or a list of vectors `list(c(min,max), c(min,max),...)` with length
#'                 `ncol(x)`, where columnwise ranges can be specified by the user.
#'              This parameter takes effect for `which='kGDM2'` and `which='kGower'`,
#'              and has no effect for `which='kModes'`.
#'              `xrange=NULL` defaults to `all` for `kGDM2`, and to `columnwise` for `kGower`.
#' @param xmethods An optional character vector of length `ncol(x)` that specifies
#'                 the distance measure for each column of `x`. Currently only used for
#'                 `kGower`. For `kGower`, `xmethods=NULL` results in the use of
#'                 default methods for each column of `x`. For more information on
#'                 allowed input values, and default measures, see the Details section.
#' @param trim See [flexclust::kccaFamily()]: Proportion of points trimmed in
#'             robust clustering.
#' @param groupFun See [flexclust::kccaFamily()]: A character string specifying
#'                 the function for clustering. Default is `'minSumClusters'`.
#'   
#' @return An object of class `flexclust::kccaFamily`. When using it within `flexclust::kcca`,
#' `kccaExtendedFamily(which='kModes')` will result in an object able to conduct
#' kModes clustering, `which='kGDM2'` will result in an object for K-centroids clustering using
#' GDM2 distance, and `which='kGower'` will result in an object for k-centroids clustering using
#' Gower's distance.
#' The output of `flexclust::kcca` will be an object of class `kcca`, and can thus
#' be used in `flexclust`'s plotting, stepwise-clustering, and bootstrapping methods.
#'    
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Gower, JC (1971).
#'   *A General Coefficient for Similarity and Some of Its Properties.*
#'   Biometrics, 27(4), 857-871.
#'   \doi{doi:10.2307/2528823}
#' - Kaufman, L, Rousseeuw, P (1990).
#'   *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics.
#'   \doi{doi:10.1002/9780470316801}
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   \doi{doi:10.1016/j.csda.2005.10.006}
#' - Walesiak, M (1993). *Statystyczna Analiza Wielowymiarowa w Badaniach Marketingowych.*
#'   Wydawnictwo Akademii Ekonomicznej, 44-46.
#' - Weihs, C, Ligges, U, Luebke, K, Raabe, N (2005).
#'   *klaR Analyzing German Business Cycles.* In: Data Analysis an
#'   Decision Support, Springer: Berlin. 335-343.
#'   \doi{doi:10.1007/3-540-28397-8_36}
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
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kGDM2'))
#' flexclust::kcca(dat, k=3, family=kccaExtendedFamily('kGDM2',
#'                                                     xrange='columnwise'))
#' # Example 3: kGower
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower'))
#' nas <- sample(c(TRUE,FALSE), prod(dim(dat)), replace=TRUE, prob=c(0.1,0.9)) |> 
#'    matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower'))
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower',
#'                                            xrange='all'))
#' flexclust::kcca(dat, 3, kccaExtendedFamily('kGower',
#'                                            xmethods=c('distEuclidean',
#'                                                       'distEuclidean',
#'                                                       'distJaccard',
#'                                                        'distManhattan',
#'                                                        'distManhattan',
#'                                                        'distSimMatch')))
#' #the case where column 2 is a binary variable, but is symmetric
#' 
#' @seealso
#' [flexclust::kcca()](https://search.r-project.org/CRAN/refmans/flexclust/html/kcca.html)
#' [flexclust::stepFlexclust()](https://search.r-project.org/CRAN/refmans/flexclust/html/stepFlexclust.html)
#' [flexclust::bootFlexclust()](https://search.r-project.org/CRAN/refmans/flexclust/html/bootFlexclust.html)
#'
#' @import flexclust

#' @export
kccaExtendedFamily <- function(which=c('kModes', 'kGDM2', 'kGower'),
                               cent=NULL, preproc=NULL,
                               xrange=NULL, xmethods=NULL,
                               trim=0, groupFun='minSumClusters') { #the last two are unused leftovers from kcca, should probably not provide them
                                 
  if(which=='kModes') {
    
    distGen <- NULL
    dstfnc <- distSimMatch
    cent <- centMode

  }
  
  if(which=='kGDM2') {
    
    if(is.null(xrange)) xrange <- 'all'
    
    rng <- .rangeMatrix(xrange)
    
    distGen <- function(x) .projectIntofx(x, rangeMatrix=rng)
    dstfnc <- distGDM2
    
    if(is.null(cent)) {
      cent <- function(x){
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
      warning('No column-wise distance measures specified, default measures
            for each column will be used.')
      distGen <- function(x) {
        #I apologize for the use of parent.frame(), but didn't know how else to fix it.
        #however I do think it's ok here because 1) 'xclass' is not a generic method
        #and 2) because I only use them once in the beginning (unlike the dists),
        #so I can't really get lost in the frames
        xcls <- get('xclass', parent.frame())
        .ChooseVarDists(xcls)
      }
      preproc <- function(x) {
        xcls <- get('xclass', parent.frame())
        .ScaleVarSpecific(x, rangeMatrix=rng,
                                    xclass=xcls)
      }
    } else {
      distGen <- function(x, ...) {
        if(!all(xmethods %in% c('distEuclidean', 'distManhattan',
                                'distSimMatch', 'distJaccard')))
          stop('Specified columnwise xmethod not implemented!')
        return(xmethods)
      }
      preproc <- function(x) .ScaleVarSpecific(x, rangeMatrix=rng,
                                               xclass=xmethods)
    }
    
    dstfnc <- distGower
    
    #the default combo in the paper for distGower was centMin. now that x is scaled in the beginning, I think centOptim is sufficient
    if(is.null(cent)) {
      cent <- function(x){
        centOptimNA(x, dist = \(y, centers) {
          distGower(y, centers, genDist=genDist)
        }) #filler cent, will be recreated in the function
      }
    }
    
  }
  
  flexclust::kccaFamily(name=which,
                        dist=dstfnc,
                        cent=cent,
                        genDist=distGen,
                        preproc=preproc,
                        trim=trim, groupFun=groupFun)
                              
}

