# Stand: 18.01.25

#old distGower_ordinal that was used in paper (written just for ordinal, not mixed)
#can be found in the AJS paper folder, is now replaced here as the other options work properly

# .ChooseVarDists: helper that maps default distances to 4 different variable types.
#                  A different option is to specifically provide a character vector
#                  of length ncol(x) to kccaFamilyGower(xmethods) that specifies the
#                  distance for each variable. In the latter case, this helper is
#                  circumvented.
#' @param xclass Character vector of length=ncol(x) of classes of variables in x,
#'               for example as obtained by sapply(data.frame(x), data.class).
#'                   Default variable specific methods will be mapped to each variable class:
#'                     - 'numeric' or 'integer': squared Euclidean distance ('distEuclidean'),
#'                                               expects data to be scaled previously (f.i. by
#'                                               .ScaleVarSpecific) 
#'                     - 'logical': Jaccard distance ('distJaccard')
#'                     - 'ordered': Manhattan distance ('distManhattan'),
#'                                               expects data to be scaled previously (f.i. by
#'                                               .ScaleVarSpecific) 
#'                     - 'factor' (i.e. categorical): Simple Matching Distance ('distSimMatch')
#'                  The treatment of the different variables can thus influenced in two ways:
#'                  First, by specifing the xmethods parameter in kccaFamilyGower(),
#'                  Secondly by their coding. F.i. for a symmetric logical variable,
#'                  symmetric treatment can be achieved by coding them either as
#'                  categorical or numeric instead of logical.
#' written after: Gower (1971), Kaufman & Rousseeuw (1990).
.ChooseVarDists <- function(xclass) {
  
  if(!is.null(dim(xclass))) {#compatibility option, where .ChooseVarDists is used outside of kcca, and directly on x
    if(is.data.frame(xclass)) {
      xclass <- sapply(xclass, data.class)
    } else {
      xclass <- rep('numeric', ncol(xclass))
    }
  }
  
  ifelse(xclass == 'numeric', 'distEuclidean', #data.class returns 'numeric' also for 'integer'
         ifelse(xclass == 'logical', 'distJaccard',
                ifelse(xclass == 'ordered', 'distManhattan',
                       ifelse(xclass == 'factor', 'distSimMatch',
                              'Error: no default method implemented for this class'
                       )
                )
         )
  ) |> setNames(names(xclass))
}


#Variable scaling after Gower and Kaufman+Rousseeuw (center around
#the minimum value, and divide by range). Only performed on numeric
#or ordinal variables. In the latter case, the function expects the
#variables to be coded numerically ranging from `1:max(level)` in steps of 1.
#' @param x a numerically coded matrix.
#' @param xclass character vector of either:
#'               - variable classes as obtained by `data.class`; or:
#'               - proposed variable-specific default distances, f.i.
#'                 as obtained by `.ChooseVarDists(x)`;
#'               with `length(xclass)==ncol(x)`
#' @param rangeMatrix expects a function that has been previously created
#'                    with `.rangeMatrix(xrange).` If it is NULL, it will be
#'                    created on `xrange`.
#' @param xrange is a compatibility parameter so the helper runs outside
#'               of the `kccaFamilyGower` concept, but within `kccaFamilyGower`,
#'               `.rangeMatrix(xrange)` is run previously
.ScaleVarSpecific <- function(x, xclass,
                              rangeMatrix=NULL, xrange=NULL) {
  
  if(is.null(rangeMatrix)) {
    rng <- .rangeMatrix(xrange)
  } else {
    rng <- rangeMatrix
  }
  
  if(is.null(colnames(x))) colnames(x) <- 1:ncol(x) #don't need this here but later in the dist #I think it's actually obsolete now
  
  cols2scl <- xclass %in% c('numeric', 'ordered', #compatibility with xclass=sapply(dat, data.class)
                            'distEuclidean', 'distManhattan') #compatibility with xmethods
  
  rng <- rng(x[, cols2scl, drop=F])
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  x[,cols2scl] <- scale(x[,cols2scl], center=rng[1,], scale=scl)
  return(x)
}

#helper function to calculate weights for Gower's distance
#' @param x a numerically coded matrix.
#' @param centers a numerically coded matrix that is compatible with,
#'                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#'                `nrow(x)>=nrow(centers)`.
#' @param distances character vector of variable specific distances,
#'                  f.i. as derived from `.ChooseVarDists(x)`. Length
#'                  needs to be equal to `ncol(x)`.
#' @return a logical array of dim `nrow(x)`X`ncol(x)`X`nrow(centers)`.
#'         `TRUE` in `ijk` if `x[i,j]` and `centers[k,j]` are both non-
#'         missing, and, in the case of logical variables calculated with
#'         `distJaccard`, not more than one variables is equal to 0; `FALSE`
#'         otherwise. 
.delta <- function(x, centers, distances) {
  K <- nrow(centers)
  delta <- sapply(1:K,
                  \(k) t(!(is.na(t(x)) | is.na(centers[k,]))),
                  simplify='array')
  distJ <- distances=='distJaccard'
  
  if(any(distJ)) {
    xJ <- x[, distJ, drop=F]
    cJ <- centers[, distJ, drop=F]
    for(k in 1:K){
      delta[, distJ, k] <- t((t(xJ) + cJ[k,])>0)
    }
    delta[which(is.na(delta))] <- FALSE
  }
  
  delta
  
}


#helper function to calculate Gower's distance on _mixed variable types_, and/or
#in the presence of missing values
#' @param x a numerically coded matrix. If `distEuclidean` or
#'          `distManhattan` are to be used on some columns, these
#'          columns need to be scaled previously, f.i. by
#'          `.ScaleVarSpecific(x)`.
#' @param centers a numerically coded matrix that is compatible with,
#'                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#'                `nrow(x)>=nrow(centers)`.
#' @param distances character vector of variable specific distances,
#'                  f.i. as derived from `.ChooseVarDists(x)`. Length
#'                  needs to be equal to `ncol(x)`.
#' @return a numeric array of dim `nrow(x)`X`ncol(x)`X`nrow(centers)`
#'         that in each `ijk` contains the distance between `x[i,j]`
#'         and `centers[k,j]`. Weighting with `delta`, and summing up
#'         over `j` happens in the next step.
#' @details NA handling: As proposed by Gower, and Kaufman+Rousseeuw. NAs
#'          are filled with placeholder values, and handled in the weighting
#'          step (Pairs of `c(x[i,j], centers[k,j])` with missing values receive
#'          weight 0, and values that are present in `c(x[i,], centers[k,])` are
#'          upweighted.
.distGower_mixedType <- function(x, centers, distances) {
  
  z <- array(0, dim=c(dim(x), nrow(centers)),
             dimnames=c(dimnames(x), list(NULL)))
  K <- nrow(centers)
  distE <- distances %in% 'distEuclidean'
  distM <- distances %in% 'distManhattan'
  distS <- distances %in% 'distSimMatch'
  distJ <- distances %in% 'distJaccard'
  if(sum(distE, distM, distS, distJ)!=ncol(x)) stop('Specified distance(s) not implemented in flexord::distGower')
  
  if(any(distE)) {
    xE <- x[, distE, drop=F]
    cE <- centers[, distE, drop=F]
    for(k in 1:K){
      z[, distE, k] <- t(sqrt((t(xE) - cE[k,])^2))
    }
  }
  if(any(distM)) {
    xM <- x[, distM, drop=F]
    cM <- centers[, distM, drop=F]
    for(k in 1:K){
      z[, distM, k] <- t(abs(t(xM)-cM[k,]))
    }
  }
  if(any(distS)) {
    xS <- x[, distS, drop=F]
    cS <- centers[, distS, drop=F]
    for(k in 1:K){
      z[, distS, k] <- t(t(xS) != cS[k,])
    }
  }
  if(any(distJ)) {
    xJ <- x[, distJ, drop=F]
    cJ <- centers[, distJ, drop=F]
    for(k in 1:K){
      z[, distJ, k] <- t((t(xJ) + cJ[k,])<2)
    }
  }

    z[which(is.na(z))] <- 1 #just a placeholder, NA cases are removed by weights==0
                            #z <- ifelse(is.na(z), 1, z)
  z
}

#helper function to calculate Gower's distance on _single variable types_,
#and no missing values.
#' @param x a numerically coded matrix. If `distEuclidean` or
#'          `distManhattan` are to be used on some columns, these
#'          columns need to be scaled previously, f.i. by
#'          `.ScaleVarSpecific(x)`.
#' @param centers a numerically coded matrix that is compatible with,
#'                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#'                `nrow(x)>=nrow(centers)`.
#' @param distances character vector of distance to be used on all variables.
#'                  Length needs to be equal to `ncol(x)`, but needs to be the
#'                  same for all columns. Can f.i. be derived from `.ChooseVarDists(x)`.
#' @return a numeric matrix of dim `nrow(x)`X`nrow(centers)`
#'         that in each `ik` contains the distance between `x[i,]`
#'         and `centers[k,]`
#' @details
#' This helper is added for the speed increase it provides when using
#' Gower's distance on all numeric or all ordered variables without
#' missing values.
#' 
#' The difference between this option, and using `distEuclidean` or 
#' `distManhattan` directly is that 1) when using these distances within
#' Gower's distance, they are scaled previously as described by Gower and
#' Kaufman+Rousseeuw, and 2) the are then divided by the number of columns
#' (as all weights will be 1). This results in a distance that ranges from 0-1.
#' 
#' Due to compatibility reasons, this function can also handle Jaccard distance,
#' #or Simple Matching Distance over all variables. However, in this case it
#' brings no advantages. Thus, it is recommended to use `flexclust::distJaccard`
#' or `distSimMatch` directly.
#' 
#' The checks for this helper to run (`length(unique(distances))==1` and
#' `!any(is.na(x))` happen outside of it.
.distGower_singleTypeNoNAs <- function(x, centers, distances) {

  dists <- unique(distances)
  p <- ncol(x)
  
  if(dists=='distEuclidean') {
    dstfnc <- flexclust::distEuclidean
  } else if(dists=='distManhattan') {
    dstfnc <- flexclust::distManhattan
  } else if(dists=='distSimMatch') {
    dstfnc <- distSimMatch
    p <- 1 #da ist colMeans schon drin
  } else if(dists=='distJaccard') {
    dstfnc <- flexclust::distJaccard
    p <- 1 #meine Vermutung ist, dass distJaccard ohne /ncol(x) für singletype noNA passt? Schließlich hat das schon einen Nenner? Und dass nur 1==1 bedacht wird, sollte da auch drin schon behandelt sein
  } else {stop('Specified distance not implemented in flexord::distGower')}

  dstfnc(x, centers)/p
}

#Function to calculate Gower's distance as written by Gower (1971) and
#Kaufman+Rousseeuw (1990), on numeric, ordered, logical, categorical,
#or mixed data sets, that possibly contain missing values.
#' @param x a numerically coded matrix. If `distEuclidean` or
#'          `distManhattan` are to be used on some columns, these
#'          columns need to be scaled previously, f.i. by
#'          `.ScaleVarSpecific(x)`.
#' @param centers a numerically coded matrix that is compatible with,
#'                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#'                `nrow(x)>=nrow(centers)`.
#' @param genDist character vector of variable specific distances to be used,
#'                  f.i. as derived from `.ChooseVarDists(x)`. Length
#'                  needs to be equal to `ncol(x)`. Can contain the options:
#'                  - `distEuclidean`: squared Euclidean distance between the
#'                                     scaled variables
#'                  - `distManhattan`: absolute distance between the scaled variables
#'                  - `distJaccard`: counts of zero if both binary variables are
#'                                   equal to 1, and 1 otherwise
#'                  - `distSimMatch`: Simple Matching Distance, i.e. the number of agreements
#'                                    between variables.
#'                These are then weighted, summed and normalized according to Gower and
#'                Kaufman+Rousseeuw.
#' @return a numeric matrix of dim `nrow(x)`X`nrow(centers)`
#'         that in each `ik` contains the distance between `x[i,]`
#'         and `centers[k,]`
#' @details
#' NA handling:
#' As proposed by Gower, and Kaufman+Rousseeuw. NAs are filled with
#' placeholder values, and handled by weighting (Pairs of
#' `c(x[i,j], centers[k,j])` with missing values receive weight 0,
#' and values that are present in `c(x[i,], centers[k,])` are upweighted.
#' 
#' This distance is designed for mixed variable types and/or data with
#' missing values. It can also bring advantages in cases of strictly numeric
#' or strictly ordered variables without missing values, but is not recommended
#' in cases of only binary or only categorical variables.
#' @export
distGower <- function(x, centers, genDist) {
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  
  if(length(unique(genDist))==1 && !any(is.na(x))) {
    
    .distGower_singleTypeNoNAs(x, centers, distances=genDist)
  
  } else {
    
    delta <- .delta(x, centers, distances=genDist)
    
    z <- .distGower_mixedType(x, centers, distances=genDist)
    
    z <- sapply(1:nrow(centers),
                \(k) rowSums(z[,,k,drop=F]*delta[,,k,drop=F])/rowSums(delta[,,k,drop=F]))
    z[which(is.nan(z))] <- 0 #aus distJaccard, dort setzt Fritz diesen Fall auch 0
    
    z
  }
}

#A wrapper for `flexclust::kccaFamily` to conduct
# K-centroids clustering with Gower's distance. It is intended for use in
# `flexclust::kcca` functions built upon it.
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
#' @return A custom `kccaFamily` object using `distGower` as the 
#'    distance function; `.ScaleVarSpecific` (scaling of numeric and ordinal
#'    variables as proposed by Gower, 1970, and Kaufman+Rousseeuw, 1990) for
#'    data preprocessing (slot `preproc`); and chooses the distances for each
#'    variable either as specified by the user in `xmethods`, or proposes default
#'    methods via `.ChooseVarDists` (slot `genDist`).
#'@export
kccaFamilyGower <- function(cent=NULL,
                            xrange='columnwise', xmethods=NULL,
                            trim=0, groupFun='minSumClusters') {
  
  rng <- .rangeMatrix(xrange)
  
#  preproc <- function(x) .ScaleGower(x, rangeMatrix=rng) #archived, simpler version, as in data wide xrange options, the scaling destroys binary variables
  
  if(is.null(xmethods)) {
    warning('No column-wise distance measures specified, default measures
            will be used. Make sure that the data object x for your clustering
            procedure is a correctly coded dataframe.') #in fact, it won't work if is.null(xmethods) && is.matrix(x)
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
  
  #the default combo in the paper for distGower was centMin. now that x is scaled in the beginning, I think centOptim is sufficient
  if(is.null(cent)) {
    cent <- function(x){
      centOptimNA(x, dist = \(y, centers) {
        distGower(y, centers, genDist=genDist)
      }) #filler cent, will be recreated in the function
    }
  }
  
  flexclust::kccaFamily(name='kGower',
                        dist=distGower,
                        genDist=distGen,
                        cent=cent,
                        preproc=preproc,
                        trim=trim, groupFun=groupFun)
}

#' @examples
#' #1) single variable type case with no missings:
#' 
#' dat <- matrix(sample(1:5, 60, replace = TRUE), nrow = 10)
#' #1.1) choose distances for each variable
#' (xcls <- .ChooseVarDists(dat)) #defaults to 'distEuclidean' for all.
#' #alternative treatment can be specified, f.i. rep('disManhattan', ncol(dat))
#' #1.2) scale according to Gower and Kaufman and Rousseeuw:
#' (dat <- .ScaleVarSpecific(dat, xclass=xcls,
#'                          xrange=c(1,6))) #let's assume that the option 6 was available for all variables, but was never used
#' #1.3) choose centers:                                     
#' initcenters <- dat[sample(1:10, 3),]
#' #1.4) calculate Gower's distance
#' distGower(dat, initcenters, genDist=xcls)
#' 
#' #calculate distance matrix (f.i. for PAM) via:
#' distGower(dat, dat, genDist=xcls) |> as.dist()
#' 
#' #within k-centroids clustering:
#' flexclust::kcca(dat, 3, kccaFamilyGower(xrange=c(1,6)))
#' #or using Manhattan distance:
#' kcca(dat, 3, kccaFamilyGower(xmethods=rep('distManhattan', 6)))
#' 
#' #2) single variable type case with missing values:
#' nas <- sample(c(T,F), prod(dim(dat)), replace=TRUE, prob=c(0.1,0.9)) |> 
#'    matrix(nrow=nrow(dat))
#' dat[nas] <- NA
#' #repeat steps as above...or just do:
#' kcca(dat, 3, kccaFamilyGower())
#' 
#' #3) mixed variable type case with no missings:
#' dat <- data.frame(cont = sample(1:100, 10, replace=T)/10,
#'                    bin_sym = as.logical(sample(0:1, 10, replace=T)),
#'                    bin_asym = as.logical(sample(0:1, 10, replace=T)),                     
#'                    ord_levmis = factor(sample(1:5, 10, replace=T),
#'                                        levels=1:6, ordered=T),
#'                    ord_levfull = factor(sample(1:4, 10, replace=T),
#'                                         levels=1:4, ordered=T),
#'                    nom = factor(sample(letters[1:4], 10, replace=T),
#'                                 levels=letters[1:4]))
#' #3.1) choose distances for each variable
#' (xcls <- .ChooseVarDists(dat))
#' #uses default distance for each variable type. If f.i. for 'bin_sym',
#' #symmetric instead of asymmetric treatment is desired, the vector needs
#' #to be created by hand (f.i. c('distEuclidean', 'distEuclidean', ...))
#' #3.2) convert to matrix, and scale:
#' datmat <- .ScaleVarSpecific(data.matrix(dat), xclass=xcls,
#'                             xrange='columnwise')
#' #caution, xrange='all' would not make sense in this case
#' #if we wanted to acknowledge that not all available levels are used in
#' #variable 'ord_levmis', we could instead provide a list of range vectors
#' #for xrange
#' #3.3) choose centers:                                     
#' initcenters <- datmat[sample(1:10, 3),]
#' #3.4) calculate Gower's distance
#' distGower(datmat, initcenters, genDist=xcls)
#' 
#' #calculate distance matrix via:
#' distGower(datmat, datmat, genDist=xcls) |> as.dist()
#' 
#' #within kcca:
#' kcca(dat, 3, kccaFamilyGower()) #using defaults for centroid, xrange, and xmethods
#'
#' #4) mixed variable type case with missing values:
#' dat[nas] <- NA
#' kcca(dat, 3, kccaFamilyGower(xrange='columnwise'))
