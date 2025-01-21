#20.01.25

#' Extending K-Centroids clustering to (mixed-with-)ordinal data
#' 
#' @description
#' This wrapper creates objects of class `flexclust::kccaFamily`,
#' which can be used within `flexclust::kcca` to conduct K-centroids
#' clustering using the following methods:
#'  - **kModes** (after Weihs et\~al., 2005)
#'  - **kGower** (Gower's distance after Kaufman \& Rousseeuw, 1990,
#'                and a specified centroid)
#'  - **kGDM2** (GDM2 distance after Walesiak et\~al., 1993, and a
#'              specified centroid)
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
#' - Numeric and/or ordinal variables are scaled by \eqn{} #hi
#' #add equation, then distance (partials+summing up), then centroid
#' #add note how colwise dists can be manipulated (after params I guess)
#' 
#' #kGDM2 
#' 
#' @param cent See [flexclust::kccaFamily()]: Function for determining cluster
#'   centroids. `NULL` (default) defaults to `flexclust::centOptim`, a general
#'   purpose optimizer.
#' @param preproc See [flexclust::kccaFamily()]: Preprocessing function to be
#' applied to the data before clustering.
#' @param trim See [flexclust::kccaFamily()]: Proportion of points trimmed in
#' robust clustering.
#' @param groupFun See [flexclust::kccaFamily()]: A character string specifying
#' the function for clustering.
#'   Default is `'minSumClusters'`.
#'   
#' @return A custom `kccaFamily` object using `distGower` as the 
#'    distance function; `.ScaleVarSpecific` (scaling of numeric and ordinal
#'    variables as proposed by Gower, 1970, and Kaufman+Rousseeuw, 1990) for
#'    data preprocessing (slot `preproc`); and chooses the distances for each
#'    variable either as specified by the user in `xmethods`, or proposes default
#'    methods via `.ChooseVarDists` (slot `genDist`).
#'    
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Kaufman, L, Rousseeuw, P (1990).
#'   *Finding Groups in Data: An Introduction to Cluster Analysis.*
#'   Wiley Series in Probability and Statistics.
#'   [doi:10.1002/9780470316801](https://doi.org/10.1002/9780470316801)
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   [doi:10.1016/j.csda.2005.10.006](https://doi.org/10.1016/j.csda.2005.10.006)
#' - Walesiak, M (1993). *Statystyczna Analiza Wielowymiarowa w Badaniach Marketingowych.*
#'   Wydawnictwo Akademii Ekonomicznej, 44-46.
#' - Weihs, C, Ligges, U, Luebke, K, Raabe, N (2005).
#'   *klaR Analyzing German Business Cycles.* In: Data Analysis an
#'   Decision Support, Springer: Berlin. 335-343.
#'   [doi:10.1007/3-540-28397-8_36](https://doi.org/10.1007/3-540-28397-8_36)


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