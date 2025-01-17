# 09.12.24

#Gower's distance

#old distGower_ordinal that was used in paper (written just for ordinal, not mixed)
#' @param xrange 'data range' for range(x); 'variable specific' for apply(x, 2, range);
#'               range vector of c(min,max); or list of range vectors for each variable
#'               (length of list == ncol(x), each item is num. vector c(min,max))
.distGower_ordinal <- function(x, centers, xrange=NULL) {
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  
  rng <- .rangeMatrix(xrange)(x)
  
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  xr <- scale(x, center=rng[1,], scale=scl)
  centr <- scale(centers, center=rng[1,], scale=scl)
  
  for(k in 1:nrow(centers)) {
    z[, k] <- colMeans(abs(t(xr) - centr[k, ]))
  }
  z  
}

#' Function can't actually be used directly on x in dist, as dist is 'data.matriced'
#' and thus loses class information. Instead, I added an if-clause in kcca, where
#' I extract the class info is extracted previously, and stored in the (at)infosOnX slot)
#' @param xclass Character vector of length=ncol(x) of EITHER:
#'                1) distances to be used for each variable. Available options are
#'                   'distEuclidean' (squared Euclidean distance), 'distManhattan'
#'                   (absolute distance), 'distJaccard' (Jaccard distance for
#'                   asymmetric binary variables), and 'distSimMatch' (simple Matching
#'                   distance, i.e. inequality between values). OR:
#'                2) classes of variables in x, for example as obtained by
#'                   sapply(data.frame(x), data.class).
#'                   Default methods will be mapped to each variable class:
#'                     - 'numeric' or 'integer': squared Euclidean distance
#'                     - 'logical': Jaccard distance
#'                     - 'ordered': Manhattan distance (after adequate preprocessing)
#'                     - 'factor' (i.e. categorical): Simple Matching Distance
#'                  The treatment of the different variables can thus also be influenced
#'                  by their coding. F.i. for a symmetric logical variable, symmetric treatment
#'                  can be achieved by coding them either as categorical or numeric instead of
#'                  logical.
#' written after: Gower (1971), Kaufman & Rousseeuw (1990).
#' Note: for by-the-book clustering with Gower's distance, scale-centering x by
#'       min and range is necessary for numeric or ordered variables. A helper
#'       for this can be found in .ScaleGower. (which applies it to all variables.
#'       BUT, scaling has no effect for logical and unordered factor variables,
#'       so no problem.)
.ChooseVarDists <- function(xclass) {
  
  if(!is.null(dim(xclass))) {#the case where the function *could* be used outside of kcca, directly on x
    #this should never evaluate to TRUE within kcca, else it'll just return
    #rep('distEuclidean', ncol(x)), and then we've got a problem
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

.rangeMatrix <- function(xrange) { #TODO: this is also created in centroidFunctions. Choose place to avoid duplicates. Either leave in centroidFunctions.R, or start the generic.R file, where I stash all .-Functions
  #xrange: response level range of the variables,
  #       implemented options:
  #       'data range': range of x, same for all variables
  #       'variable specific': range of each variable in the data
  #       c(lower, upper): range vector specified by user, upper and lower limit for all variables
  #       list(x1= c(lower, upper), x2=c(lower, upper), ...): list of range vectors specified by user, upper and lower limits are variable specific
  if(all(xrange=='data range')) {
    rng <- function(x) {
      rep(range(x, na.rm=T), ncol(x)) |>
        matrix(nrow=2)
    }
  } else if(all(xrange=='variable specific')) {
    rng <- function(x) {
      apply(x, 2, range, na.rm=T)
    }
  } else if(is.vector(xrange, mode='numeric')) {
    if(length(xrange) != 2)
      stop('Either supply 1 range vector, or list of ranges for all variables')
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
#'                    with .rangeMatrix(xrange). If it is NULL, see the next param.
#' @param xrange is a compatibility parameter so the helper runs outside
#'               of the kccaFamilyGower concept, but within kccaFamilyGower,
#'               .rangeMatrix(xrange) is run previously
.ScaleVarSpecific <- function(x, xclass,
                              rangeMatrix=NULL, xrange=NULL) {
  
  if(is.null(rangeMatrix)) {
    rng <- .rangeMatrix(xrange)
  } else {
    rng <- rangeMatrix
  }
  
  if(is.null(colnames(x))) colnames(x) <- 1:ncol(x) #don't need this here but later in the dist
  
  cols2scl <- xclass %in% c('numeric', 'ordered', #compatibility with xclass=sapply(dat, data.class)
                            'distEuclidean', 'distManhattan') #compatibility with xmethods
  
  rng <- rng(x[, cols2scl, drop=F])
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  x[,cols2scl] <- scale(x[,cols2scl], center=rng[1,], scale=scl)
  return(x)
}

##this is the option that just scales over all, cause it doesn't matter for binary and nominal:
#NO, ERROR! It doesn't matter if range is calculated variable specific,
#or if user specifies for each variable. HOWEVER, it does destroy distJaccard
#for ranges that cover the whole data set, variables aren't 0,1 anymore afterward
#NEED TO REACTIVATE VARIABLE SPECIFIC SCALING, this is archived
#' @param rangeMatrix expects a function that has been previously created
#'                    with .rangeMatrix(xrange). If it is NULL, see the next param.
#' @param xrange is a compatibility parameter so the helper runs outside
#'               of the kccaFamilyGower concept, but within kccaFamilyGower,
#'               .rangeMatrix(xrange) is run previously
.ScaleGower <- function(x, rangeMatrix=NULL,
                        xrange=NULL) {

  if(is.null(rangeMatrix)) {
    rng <- .rangeMatrix(xrange)
  } else {
    rng <- rangeMatrix
  }
  
  if(is.null(colnames(x))) colnames(x) <- 1:ncol(x) #don't need this here but later in the dist
  
  rng <- rng(x)
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  scale(x, center=rng[1,], scale=scl)
}

#helper function to calculate Gower's weights
#     (i.e. both values non-missing, and in the case of distJaccard,
#      not both values =0)
#' @param distances expects genDist object, f.i. as derived from .ChooseVarDists(x)
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


#helper function to Gower's distance on _mixed vartypes_
#     (i.e. x needs to be scaled previously, distance function is
#      chosen for each column, NAs are handled as: omitted, nonmissing vars are upweighted).
#      output: array of dim nrow(x) X ncol(x) X nrow(centers)
#               (Weighting with delta, and summing up happens at the end of the distGower function)
#' @param distances expects genDist object, f.i. as derived from .ChooseVarDists(x)
.distGower_mixedType <- function(x, centers, distances) {
  
  #  dists <- unique(distances)
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

#helper function to Gower's distance on a _single vartype_ with no NAs
#checks need to happen outside of it
#     (i.e. x needs to be scaled previously, distance function is
#      chosen for the single vartype).
#      output: array of dim nrow(x) X ncol(x) X nrow(centers)
#               (Weighting with delta, and summing up happens at the end of the distGower function)
#' @param distances expects genDist object, f.i. as derived from .ChooseVarDists(x)
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


#' @param genDist i.e. genDist expects char.vector with dist to be used
#'                 for each variable, as for example created by .ChooseVarDists
#' @param x, centers: dist expects these to be already scaled 'gowerly'
distGower <- function(x, centers, genDist) {
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  
  #if single vartype, and no NAs --> z=simpledistancemethod/ncol(x)
  #for distEuclidean, distManhattan (as in our study), this still follows
  #Gower, as the variables have been scaled previously, and are didided by
  #the number of parameters (i.e. dist will always range between 0 and 1).
  #For distJaccard and distSimMatch it is redundant, but I still included
  #them for now
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




kccaFamilyGower <- function(cent=NULL,
                            xrange=NULL, xmethods=NULL,
                            trim=0, groupFun='minSumClusters') {
  
  rng <- .rangeMatrix(xrange)
  
  preproc <- function(x) .ScaleGower(x, rangeMatrix=rng)
  #if I were to use .ScaleVarSpecific, this would get
  #more complicated cause in case of is.null(xmethods),
  #I'd need to access the primed family to get xclass
  
  
  if(is.null(xmethods)) {
    warning('No column-wise distance measures specified, default measures
            will be used. Make sure that the data object x for your clustering
            procedure is a correctly coded dataframe.') #in fact, it won't work if is.null(xmethods) && is.matrix(x)
    distGen <- function(x, ...) { #TODO: remove ... when GDM2's genDist is fixed
#      if('xclass' %in% names(xclass)) xclass <- xclass$xclass #again, the necessary catcher for the family@infosOnX
#      .ChooseVarDists(xclass=xclass) #adding x here cause that's the setup I need for the relevant code line, but fun is independent of x (at this stage, x is already 'numericized')
      #now: trying to create a function factory, so that when calling family@genDist in kcca, it'll just create the function (x doesn't do anything), and the function then has access to the xclass object that was created a few lines ago in kcca
      #not successfull, had to work with parent.frame(). I think it's okay in this case though, as this function is only used in this instance, and not other, lower-lying ones. Also, it's not supposed to run outside of kcca, right? (if one wanted to use it, could just use .ChooseVarDists)
      xcls <- get('xclass', parent.frame())
      .ChooseVarDists(xcls)
    }
  } else {
    distGen <- function(x, ...) { #TODO: dots necessary right now, remove when distGDM2 is updated
#      if('xmethods' %in% names(xmethods)) xmethods <- xmethods$xmethods #s.o.
      if(!all(xmethods %in% c('distEuclidean', 'distManhattan',
                              'distSimMatch', 'distJaccard')))
        stop('Specified columnwise xmethod not implemented!')
      return(xmethods)
    }
  }
  
  # if(is.null(cent)) {
  #   cent <- function(x) {
  #     #if(is.function(xmethods)) xmethods <- xmethods(x) #this is now handled directly in dist
  #     centMin(x, xrange=xrange,
  #             dist = \(y, centers) {
  #               distGower(y, centers,
  #                         genDist=xmethods)
  #     })
  #   }
  # }
  #now that the variables are scaled, is centOptim sufficient?
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
                        xrange=rng, xmethods=xmethods,
                        trim=trim, groupFun=groupFun)
}

kccaFamilyGowerTestDataRange <- function(cent=NULL,
                            xrange=NULL, xmethods=NULL,
                            trim=0, groupFun='minSumClusters') {
  
  rng <- .rangeMatrix(xrange)
  
#  preproc <- function(x) .ScaleGower(x, rangeMatrix=rng)
  #archived, simpler version, as in the case of xrange='data range' or xrange=c(min,max)
  #the scaling destroys binary variables
  
  
  if(is.null(xmethods)) {
    warning('No column-wise distance measures specified, default measures
            will be used. Make sure that the data object x for your clustering
            procedure is a correctly coded dataframe.') #in fact, it won't work if is.null(xmethods) && is.matrix(x)
    distGen <- function(x, ...) { #TODO: remove ... when GDM2's genDist is fixed
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
    distGen <- function(x, ...) { #TODO: dots necessary right now, remove when distGDM2 is updated
      #      if('xmethods' %in% names(xmethods)) xmethods <- xmethods$xmethods #s.o.
      if(!all(xmethods %in% c('distEuclidean', 'distManhattan',
                              'distSimMatch', 'distJaccard')))
        stop('Specified columnwise xmethod not implemented!')
      return(xmethods)
    }
    preproc <- function(x) .ScaleVarSpecific(x, rangeMatrix=rng,
                                             xclass=xmethods)
  }
  
  #the default combo in the paper for distGower was centMin.
  #now that x is scaled in the beginning, I think centOptim is
  #sufficient
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
                        xrange=rng, xmethods=xmethods,
                        trim=trim, groupFun=groupFun)
}

if(FALSE){ #Example zone
  #mixed data case
  (dat <- data.frame(cont = sample(1:100, 10, replace=T)/10,
                     bin_sym = as.logical(sample(0:1, 10, replace=T)),
                     bin_asym = as.logical(sample(0:1, 10, replace=T)),                     
                     ord_levmis = factor(sample(1:5, 10, replace=T),
                                         levels=1:6, ordered=T),
                     ord_levfull = factor(sample(1:4, 10, replace=T),
                                          levels=1:4, ordered=T),
                     nom = factor(sample(letters[1:4], 10, replace=T),
                                  levels=letters[1:4])))
  
  datrng <- apply(datmat, 2, range) |> data.frame() |> 
    as.list()
  datmat <- data.matrix(dat) |> 
    .ScaleGower(xrange='variable specific')
  datcent <- datmat[sample(1:nrow(dat), 3, replace=F),]
  
  #mixed data with NAs
  datNA <- dat
  nas <- sample(c(T,F), prod(dim(dat)),
                replace=T, prob=c(0.1,0.9))
  datNA[nas] <- NA
  datNAmat <- data.matrix(datNA) |> 
    .ScaleGower(xrange=datrng)
  datcentNA <- datmatNA[sample(1:nrow(dat), 3, replace=F),]

  #hi
  #(hopefully) correct .ScaleVarSpecific exists now,
  #test it, and test clustering with kccaFamilyGowerTestDataRange
  
  genDist <- .ChooseVarDists(dat)
  
  datmatNAscld <- .ScaleGower(datmatNAscld, xrange=datrng)
  datcentNAscld <- .ScaleGower(datcentNA, xrange=datrng)
  
  #single variable type case
  data('risk', package='MSA')
  riskcent <- risk[sample(1:nrow(risk), 3, replace=F),]
  
  
  datmatscld |> distGower(datcent)

  datmatrixscaled |> distGower(centnoNAscaled, .ChooseVarDists(dat))
  datmatrixscaled |> distGower(datmatrix, .ChooseVarDists(dat)) |> as.dist()
  #testing the 1-vartype case (with Euclidean distance)
  xscld <- .ScaleGower(x)
  centscld <- .ScaleGower(centers)
  distGower(xscld, centscld, .ChooseVarDists(x))
  distGower(xscld, xscld, .ChooseVarDists(x)) |> as.dist()
  #testing the 1-vartype case (with Jaccard)
  distGower(xscld, centscld, rep('distJaccard', 6))
  #testing the nrow(centers)=1 case
  distGower(datmatrixscaled, centnoNAscaled[1,,drop=F], .ChooseVarDists(dat))
  #testing the nrow(x) && nrow(centers) == 1 case
  (h <- datmatrixscaled[1,,drop=F])
  distGower(h, h, .ChooseVarDists(dat))
  distGower(x[1,,drop=F], x[1,,drop=F], .ChooseVarDists(x))
  
  #testing the NA case
  distGower(datNAmatrixscaled, centNAscaled, .ChooseVarDists(dat))
  #-within kcca, stepFlexclust, bootFlexclust
  # *for the cases: mixed vars, one vartype, with centOptim (1 var case usually
  # happens for the datmatrix object), with NAs
  kcca(datmatrix, k, family=kccaFamilyGower()) #worked for now
  stepFlexclust(datmatrix, k=2:4, family=kccaFamilyGower())
  bootFlexclust(datmatrix, k=2:4, nboot=3,
                family=kccaFamilyGower())
}
