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
#'                   sapply(as.data.frame(x), \(y) paste(class(y), collapse=' ')).
#'                   Default methods will be mapped to each variable class:
#'                     - 'numeric' or 'integer': squared Euclidean distance
#'                     - 'logical': Jaccard distance
#'                     - 'ordered factor': Manhattan distance (after adequate preprocessing)
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
    if(!is.data.frame(xclass)) {
      xclass <- rep('numeric', ncol(xclass))
    } else {
      xclass <- sapply(xclass, \(y) paste(class(y), collapse=' '))
    }
  }
  
  ifelse(xclass %in% c('numeric', 'integer'), 'distEuclidean',
         ifelse(xclass == 'logical', 'distJaccard',
                ifelse(xclass == 'ordered factor', 'distManhattan',
                       ifelse(xclass == 'factor', 'distSimMatch',
                              'Error: no default method implemented for this class'
                       )
                )
         )
  ) |> setNames(names(xclass))
}

.rangeMatrix <- function(xrange) {
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

.ScaleVarSpecific <- function(x, xrange, xmethods) {
  #making scale Y/N dependent on vartype. BUT I could just scale over all
  #vartypes, it wouldn't change the resulting distances.
  #x...data matrix
  #xrange...same options as before, will be handled by .rangeMatrix
  #xmethods...character vector of length=ncol(x). Either dists to be
  #           applied to each variable, or variable class (latter one
  #           will be converted to char. vector of dists by .ChooseVarDists)
  if(!all(startsWith(xmethods, 'dist'))) xmethods <- .ChooseVarDists(xclass=xmethods)
 
  cols2scale <- xmethods %in% c('distEuclidean', 'distManhattan')
  
  rng <- .rangeMatrix(xrange)(x[,cols2scale])
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  x[,cols2scale] <- scale(x[,cols2scale], center=rng[1,], scale=scl)
  return(x)
}

.ScaleGower <- function(x, rangeMatrix=NULL,
                        xrange=NULL) {
  #this is the option that just scales over all, cause it doesn't matter for binary and nominal:
  #' @param rangeMatrix expects a function that has been previously created
  #'                    with .rangeMatrix(xrange). If it doesn't exist, see the next param.
  #' @param xrange is a compatibility parameter so the helper runs outside
  #'               of the kccaFamilyGower concept, but within kccaFamilyGower,
  #'               .rangeMatrix(xrange) will be run previously
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

#' @param genDist i.e. genDist expects char.vector with dist to be used
#'                 for each variable, as for example created by .ChooseVarDists
#' @param x, centers: dist expects these to be already scaled 'gowerly'

distGower <- function(x, centers, genDist) {
  
  #would it be sufficient here to just call flexclust::distEuclidean etc.
  #on the relevant columns, replace NAs, and sum up the different parts?
  #no: too many NAs in the end (after all, not only 1 spot NA, but the
  #whole row xdistEuclidean)
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  
  xsep <- sapply(unique(genDist), \(y) {
    x[,which(genDist==y), drop=F]
  }, simplify = F)
  centsep <- sapply(unique(genDist), \(y) {
    centers[,which(genDist==y), drop=F]
  }, simplify = F)

  #calculating the weights (used both in numer and denom)  
  delta <- sapply(1:nrow(centers),
                  \(y) !(is.na(x) | is.na(y)),
                  simplify = 'array')
  #weight adjustment for Jaccard case happens below
  
  hlp <- xsep
  
  for(dist_type in unique(genDist)) {
    
    dist_cols <- sum(genDist==dist_type)
    
    if(dist_type == 'distEuclidean') {
      hlp[[dist_type]] <- sapply(1:nrow(centers), \(i) { #leaving this '1:nrow(centers)' instead of just centsep$distManhatttan, cause I guess that makes it more compatible to the others #eh, it's all the same
        d <- sqrt((xsep[[dist_type]] - matrix(centsep[[dist_type]][i,],
                                              nrow=nrow(x),
                                              ncol=dist_cols,
                                              byrow=T))^2)
        ifelse(is.na(d), 1, d) #using 1 (=max.dist for the scaled vars) as placeholder
      }, simplify = 'array')
    }
    
    if(dist_type == 'distManhattan') {
      hlp[[dist_type]] <- sapply(1:nrow(centers), \(i) {
        d <- abs(xsep[[dist_type]] - matrix(centsep[[dist_type]][i,],
                                             nrow=nrow(x),
                                             ncol=dist_cols,
                                             byrow=T))
        ifelse(is.na(d), 1, d) #s.o.
      }, simplify = 'array')
    }
  
    if(dist_type == 'distSimMatch') {
      hlp[[dist_type]] <- sapply(1:nrow(centers), \(i) {
        d <- xsep[[dist_type]] != matrix(centsep[[dist_type]][i,],
                                         nrow=nrow(x),
                                         ncol=dist_cols,
                                         byrow=T)
        ifelse(is.na(d), 1, d)
      }, simplify = 'array')
    }
    
    if(dist_type == 'distJaccard') {
      
      centrep <- sapply(1:nrow(centers), \(i) {
        matrix(centsep[[dist_type]][i,],
               nrow=nrow(x),
               ncol=dist_cols,
               byrow=T)
      }, simplify = 'array')
      
      hlp[[dist_type]] <- sapply(1:nrow(centers), \(i) {
        d <- xsep[[dist_type]] != centrep[,,i]
        d <- ifelse(is.na(d), 1, d) #das sollte im nächsten Schritt eh auch
        #mitgehandelt werden, aber zur Sicherheit explizit
        d[xsep[[dist_type]]==0] <- 1 #außerdem könnte ich mir den Schritt eh
        #sparen, die Info muss ja sowieso ins delta
        d
      }, simplify = 'array')
      
      delta <- sapply(1:nrow(centers), \(i) {
        dlt <- xsep[[dist_type]]+centrep[,,i]
        dlt <- ifelse(is.na(dlt), 0, dlt)
        delta[,genDist==dist_type, i][dlt==0] <- 0
        delta[,,i]
      }, simplify = 'array')
    }
  }
  
  hlp <- sapply(1:nrow(centers), \(i) {
    arr <- lapply(hlp, \(y) y[,,i, drop=F])
    dims <- sapply(lapply(arr, dim), `[`, 2)
    
    arr[dims==1] <- lapply(arr[dims==1],
                           \(y) matrix(y, dimnames=dimnames(y)[1:2]))
    arr[dims>1] <- lapply(arr[dims>1], \(y) y[,,1])
    do.call(cbind, arr)[,colnames(x)]
  }, simplify = 'array')

  z <- sapply(1:nrow(centers),
              \(i) rowSums(hlp[,,i]*delta[,,i])/rowSums(delta[,,i]))
  return(z)
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
    distGen <- function(x, xclass) {
      if('xclass' %in% names(xclass)) xclass <- xclass$xclass #again, the necessary catcher for the family@infosOnX
      .ChooseVarDists(xclass=xclass) #adding x here cause that's the setup I need for the relevant code line, but fun is independent of x (at this stage, x is already 'numericized')
    }
  } else {
    distGen <- function(x, xmethods) {
      if('xmethods' %in% names(xmethods)) xmethods <- xmethods$xmethods #s.o.
      if(!all(xmethods %in% c('distEuclidean', 'distManhattan',
                              'distSimMatch', 'distJaccard')))
        stop('Specified distance not implemented in xmethod!')
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

if(FALSE){
  (dat <- data.frame(cont = sample(1:100, 10, replace=T)/10,
                     bin_sym = as.logical(sample(0:1, 10, replace=T)),
                     bin_asym = as.logical(sample(0:1, 10, replace=T)),                     
                     ord_levmis = factor(sample(1:5, 10, replace=T),
                                         levels=1:6, ordered=T),
                     ord_levfull = factor(sample(1:4, 10, replace=T),
                                          levels=1:4, ordered=T),
                     nom = factor(sample(letters[1:4], 10, replace=T),
                                  levels=letters[1:4])))
  datmatrix <- data.matrix(dat)
  
  .ScaleVarSpecific(datmatrix,
                    'variable specific',
                    xmethods=sapply(dat, \(y) paste(class(y), collapse=' ')))
  #trying out overriding colwise distances
  .ScaleVarSpecific(datmatrix,
                    'variable specific',
                    xmethods=rep('distManhattan', ncol(datmatrix)))
  #only column 'nom' looks different now. However, this should not
  #matter for the analysis because, as stated identical(1,1) gives
  #the same as identical(1/x,1/x).
  #let's see if there are runtime differences:
  microbenchmark::microbenchmark(
    colwise=.ScaleVarSpecific(datmatrix,
                      'variable specific',
                      xmethods=sapply(dat, \(y) paste(class(y), collapse=' '))),
    overall=.ScaleVarSpecific(datmatrix,
                      'variable specific',
                      xmethods=rep('distManhattan', ncol(datmatrix)))
  ) #overall clearly faster, but what about bigger objects?
  
  xmethods <- sapply(dat, \(y) paste(class(y), collapse = ' '))
  system.time(
    .ScaleVarSpecific(datmatrix[rep(1:nrow(dat), 10000),],
                      'variable specific',
                      xmethods=xmethods)
  )
  system.time(
    .ScaleGower(datmatrix[rep(1:nrow(dat), 10000),],
                'variable specific')
  ) #there .ScaleGower is slower...
  #what is critical value?
  system.time(
    .ScaleVarSpecific(datmatrix[rep(1:nrow(dat), 1000),],
                      'variable specific',
                      xmethods=xmethods)
  ); system.time(
    .ScaleGower(datmatrix[rep(1:nrow(dat), 1000),],
                'variable specific')
  ) #there is a *tiny* difference already...
  #I don't know...
  
  #I guess I'll leave the helper .ScaleVarSpecific in, and use
  #.ScaleGower for now, cause it's *considerably* easier
  
  genDist <- .ChooseVarDists(dat)
  
  #add some NAs to dat
  d <- prod(dim(dat))
  na <- vector('logical', length=d) |> 
    matrix(nrow=nrow(dat))
  na[sample(1:d, d/10, replace=F)] <- T
  datNA <- dat
  datNA[na] <- NA
  datNAmatrix <- data.matrix(datNA)
  centNA <- datNAmatrix[sample(1:nrow(dat), 3, replace=F),]
  
  datNAmatrixscaled <- .ScaleGower(datNAmatrix, xrange='variable specific')
  centNAscaled <- datNAmatrixscaled[sample(1:nrow(dat), 3, replace=F),]
  distGower(datNAmatrixscaled, centNAscaled,
            genDist=.ChooseVarDists(datNA))
  #ERMAGAAAAAAAAWD IT WORKS!!!!!!!!!!!!!!!!!
  #now I still need to clean up the accesses within kcca
  #but that's for tomorrow, my brain is dead
  distGower(x, centers, genDist=.ChooseVarDists(x))
  #ofc no workey, still need to change it so that it'll only run when 
  #the dist is used within genDist
  distGowernew(datNAmatrixscaled, centNAscaled,
               genDist=.ChooseVarDists(datNA))
  distGowernew(.ScaleGower(x),
               .ScaleGower(centers), genDist=.ChooseVarDists(x))
  #yas!
  #now replacing distGower with distGowernew
  #checken, ob das das gleiche Ergebnis bringt wie die anderen
  distEuclidean(.ScaleGower(x), .ScaleGower(centers))
  #noup, these are in fact different results and get bigger than 1
  identical(distGower(.ScaleGower(x), .ScaleGower(centers),
               genDist=rep('distManhattan', ncol(x))),
            .distGower_ordinal(x, centers))
  #yas!
  microbenchmark::microbenchmark(
    distGower(.ScaleGower(x), .ScaleGower(centers),
                 genDist=rep('distManhattan', ncol(x))),
    .distGower_ordinal(x, centers)
  ) #uiuiuiuiui das ist schon deeeutlich langsamer
  
  kcca(x, k, family=kccaFamilyGower(xmethods=rep('distManhattan', ncol(x))))
  kcca(x, k, family=kccaFamilyGower())
  kcca(dat, k, family=kccaFamilyGower())
  kcca(datNA, k, family=kccaFamilyGower())
  
  #Error occurs within the allcent step
  #family@allcent at this step is:
  function(x, cluster, k=max(cluster, na.rm=TRUE))
  {
    centers <- matrix(NA, nrow=k, ncol=ncol(x))
    for(n in 1:k){
      if(sum(cluster==n, na.rm=TRUE)>0){
        centers[n,] <- z@cent(x[cluster==n,,drop=FALSE])
      }
    }
    centers
  }
  #The error is:
  #get('z', environment(family@allcent))@cent(x[cluster==1,,drop=F])
  #Error in y[, , i, drop = F] : incorrect number of dimensions
  #and z@cent at this step is:
  #get('z', environment(family@allcent))@cent
  function(x) {
    eval(bquote({
      .(origCent)
    }))
  }
 # <environment: 0x563599f12d20>
  #and this environment contains:
  #ls(environment(get('z', environment(family@allcent))@cent))
  #[1] "centers"   "cluster"   "clustold"  "control"   "distmat"   "family"   
  #[7] "genDist"   "group"     "iter"      "k"         "MYCALL"    "N"        
  #[13] "origCent"  "origDist"  "sannprob"  "save.data" "simple"    "weights"  
  #[19] "x" 
  #with origCent being:
  #get('origCent', environment(get('z', environment(family@allcent))@cent))
  {
    centOptimNA(x, dist = function(y, centers) {
      distGower(y, centers, genDist = genDist)
    })
  }
  #and genDist being:
  #get('genDist', environment(get('z', environment(family@allcent))@cent))
  #cont         bin_sym        bin_asym      ord_levmis     ord_levfull 
  #"distEuclidean"   "distJaccard"   "distJaccard" "distManhattan" "distManhattan" 
  #nom 
  #"distSimMatch"
  
  #so let's do:
  centOptimNA(datmatrix, dist= \(y, centers) {
    distGower(y, centers, genDist = c('distEuclidean',
                                      'distJaccard',
                                      'distJaccard',
                                      'distManhattan',
                                      'distManhattan',
                                      'distSimMatch'))
  })
  
}
