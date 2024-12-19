# 06.12.24

#' Helper functions for K-centroids clustering using GDM2 distance
#'
#'@usage
#' distGDM2(x, centers, genDist=NULL, xrange=NULL)
#' .projectIntofx(x, xrange=NULL)
#' kccaFamilyGDM2(cent=NULL, preproc=NULL,
#'                xrange=NULL, xmethods=NULL,
#'                trim=0, groupFun='minSumClusters')
#'
#'@description
#' `distGDM2` implements the GDM2 distance first presented by Walesiak (1993)
#' and extended to usage in K-centroids clustering in Ernst et.~al (2025).
#' 
#' `kccafamilyGDM2` is a wrapper for `flexclust::kccaFamily` to conduct
#' K-centroids clustering with GDM2 distance. It is intended for use in
#' `flexclust::kcca` functions built upon it.
#' 
#' @details
#' `.projectIntofx` is a helper function factory that creates a function that
#' will project an arbitrary object into the space of the empirical PDF and CDF
#' of a matrix x.
#' 
#' If `distGDM2` is used within `flexclust:kcca`, the `genDist` function
#' (such as `.projectIntofx`) will be automatically primed on `x` within
#' `kcca` previously to applying the distance function.
#' If `distGDM2` is to be used on its own outside of kcca, priming will
#' occur within `distGDM2`.
#'
#' @param x A matrix of numerically coded ordinal data points. The ordinal
#' variables need to be coded as `1:length(levels(x[,i]))`.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#' `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' Within `kcca`, this is filled with the cluster centers.
#' @param genDist Function for creating a distance function that will be
#'    primed on `x`, such as f.i. `.projectIntofx`. For more information, see
#'    the 'Details' section.
#' @param xrange Range specification for the variables. Possible values are:
#'   - `NULL` (default): defaults to `'data range'`.
#'   - `'data range'`: Uses the range of the entire dataset.
#'   - `'variable specific'`: Uses column-specific ranges.
#'   - A numeric vector: Applies the same range to all columns.
#'   - A list of numeric vectors: Specifies column-specific ranges.
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
#' @return
#' - `distGDM2`: A distance matrix for each row in `x` from each row in `centers`
#'    with dimensions `c(nrow(x), nrow(centers)`.
#' - `.projectIntofx`: A function with parameter `new_x`that will project the
#'    new data object `new_x` into the space of `x`'s epdf, ecdf, and \( \tilde{F} \).
#' - `kccaFamilyGDM2`: A custom `kccaFamily` object using `distGDM2` as the 
#'    distance function, and `.projectIntofx` as the belonging distance generating
#'    function. To be used within `flexclust::kcca` or functions that build upon it.
#'
#' @examples
#' ## Example usage
#' #creating a distance matrix between two matrices based on GDM2 distance
#' dat <- matrix(sample(1:5, 50, replace = TRUE),
#'               nrow = 10, ncol = 5)
#' initcenters <- dat[sample(1:10, 3),]
#' distGDM2(dat, initcenters, genDist=.projectIntofx)
#' 
#' #a classical distance matrix for a single matrix can be obtained via:
#' as.dist(distGDM2(dat, dat, genDist=.projectIntofx))
#' 
#' #K-centroids clustering using GDM2 distance
#' flexclust::kcca(dat, k=3, family=kccaFamilyGDM2())
#' 
#' @seealso - [flexclust::kcca()](https://cran.r-project.org/package=flexclust)
#'
#' @references
#' - Leisch, F (2006). *A Toolbox for K-Centroids Cluster Analysis.*
#'   Computational Statistics and Data Analysis, 17(3), 526-544.
#'   [doi:10.1016/j.csda.2005.10.006](https://doi.org/10.1016/j.csda.2005.10.006)
#' - Walesiak, M (1993). *Statystyczna Analiza Wielowymiarowa w Badaniach Marketingowych.*
#'   Wydawnictwo Akademii Ekonomicznej, 44-46.
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#'
#' @export

distGDM2 <- function(x, centers, genDist, xrange=NULL) {
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  
  if(is.function(genDist)) {
    if(length(formals(genDist))>1) {#i.e. if genDist is 'unprimed', i.e. the function runs outside of the kcca environment
      genDist <- genDist(x, xrange)
    }
  } else { stop('No valid genDist object found') }
  
  N <- nrow(x)
  fx <- get('fx', envir=environment(genDist))
  fc <- genDist(centers)
  
  for(k in 1:nrow(centers)) {
    for(i in 1:N) {
      deltaeq <- x[i,]==centers[k,]
      numneq <- (!deltaeq)*(1/N + 2*abs(fx$Ftildex[i,]-fc$Ftildex[k,]))
      numeq <- deltaeq*fx$epdf[i,]
      denom <- sqrt(sum(1-fx$epdf[i,])*sum(1-fc$epdf[k,]))
      
      z[i,k] <- (1-(sum(1-numneq-numeq)/denom))/2
    }
  }
  z 
}

#' @export
.projectIntofx <- function(x, xrange=NULL){
  if(is.null(xrange)) xrange <- 'data range'
  
  if('xrange' %in% names(xrange)) xrange <- xrange$xrange
  #this is for the case within kcca where I pass genDist(x, family@infosOnX)
  #to make that step there more generalizable, so I can use the call in kcca
  #for different functions with different formals
  
  rng <- .rangeMatrix(xrange=xrange)(x) #code for this lies in centroids.R
  
  hats <- lapply(1:ncol(x), function(y) {
    level <- factor(x[,y], levels=seq(rng[1,y], rng[2,y]))
    pdf <- table(level)/nrow(x)
    pdf <- as.data.frame.table(pdf)
    pdf$level <- as.numeric(pdf$level)
    
    epdf <- function(i) {
      ind <- which(pdf[,'level'] <= i)
      if(length(ind) == 0) {
        return(0)
      } else {
        return(pdf$Freq[ind[length(ind)]])
      }
    }
    list(epdf=epdf, ecdf=ecdf(x[,y]))
  })
  
  names(hats) <- colnames(x)
  
  projectNew <- function(new_x) {
    fnew_x <- sapply(c('epdf', 'ecdf'), function(type) {
      z <- matrix(0, nrow=nrow(new_x), ncol=ncol(new_x),
                  dimnames=dimnames(new_x))
      for(j in 1:ncol(new_x)) {
        for(i in 1:nrow(new_x)) {
          z[i,j] <- hats[[j]][[type]](new_x[i,j])
        }
      }
      z
    }, simplify = FALSE)
    fnew_x$Ftildex <- fnew_x$ecdf - (fnew_x$epdf/2)
    return(fnew_x)
  }
  
  fx <- projectNew(x)
  
  return(projectNew)
}

#' @import flexclust
#' @export
kccaFamilyGDM2 <- function(cent=NULL, preproc=NULL,
                           xrange=NULL, xmethods=NULL,
                           trim=0, groupFun='minSumClusters') {
  if(is.null(cent)) {
    cent <- function(x){
    #genDist <- .projectIntofx(x, xrange=xrange)
    flexclust::centOptim(x, dist = \(y, centers) {
      distGDM2(y, centers, genDist=genDist)
      }) #filler cent, will be recreated in the function
    }
  }
  flexclust::kccaFamily(name='kGDM2',
                        dist=distGDM2,
                        genDist=.projectIntofx,
                        cent=cent, preproc=preproc,
                        xrange=xrange, xmethods=xmethods,
                        trim=trim, groupFun=groupFun)
}
