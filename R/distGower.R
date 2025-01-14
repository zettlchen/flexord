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
  
  k <- nrow(centers)
  typs <- unique(genDist)
  
  if(length(typs)==1) { #the single vartype case
    if(all(typs=='distEuclidean')) {
      distnc <- function(x, centers) flexclust::distEuclidean(x, centers)
    } else if(all(typs=='distManhattan')) {
      distnc <- function(x, centers) flexclust::distManhattan(x, centers)
    } else if(all(typs=='distJaccard')) {
      distnc <- function(x, centers) flexclust::distJaccard(x, centers)
    } else if(all(typs=='distSimMatch')) {
      distnc <- function(x, centers) distSimMatch(x, centers)
    } else {stop('Specified distance not implemented in flexord::distGower')}
    
    distnc(x, centers)
    
  } else {
    xsep <- sapply(typs, \(y) {
      x[,which(genDist==y), drop=F]
    }, simplify = F)
    centsep <- sapply(typs, \(y) {
      centers[,which(genDist==y), drop=F]
    }, simplify = F)
    
    #calculating the weights (used both in numer and denom)  
    delta <- sapply(1:k,
                    \(y) !(is.na(x) | is.na(y)),
                    simplify = 'array')
    #weight adjustment for Jaccard case happens below
    
    hlp <- xsep
    
    for(dist_type in typs) {
      
      dist_cols <- sum(genDist==dist_type)
      
      if(dist_type == 'distEuclidean') {
        hlp[[dist_type]] <- sapply(1:k, \(i) { 
          d <- sqrt((xsep[[dist_type]] - matrix(centsep[[dist_type]][i,],
                                                nrow=nrow(x),
                                                ncol=dist_cols,
                                                byrow=T))^2)
          ifelse(is.na(d), 1, d) #using 1 (=max.dist for the scaled vars) as placeholder
        }, simplify = 'array')
        if(k == 1) {
          hlp[[dist_type]] <- array(hlp[[dist_type]],
                                    dim=c(nrow(x), dist_cols, k),
                                    dimnames=list(NULL,
                                                  names(which(dist_type==genDist)), NULL))
        }
      }
      
      if(dist_type == 'distManhattan') {
        hlp[[dist_type]] <- sapply(1:k, \(i) {
          d <- abs(xsep[[dist_type]] - matrix(centsep[[dist_type]][i,],
                                              nrow=nrow(x),
                                              ncol=dist_cols,
                                              byrow=T))
          ifelse(is.na(d), 1, d) #s.o.
        }, simplify = 'array')
        if(k == 1) {
          hlp[[dist_type]] <- array(hlp[[dist_type]],
                                    dim=c(nrow(x), dist_cols, k),
                                    dimnames=list(NULL,
                                                  names(which(dist_type==genDist)), NULL))
        }
      }
      
      if(dist_type == 'distSimMatch') {
        hlp[[dist_type]] <- sapply(1:k, \(i) {
          d <- xsep[[dist_type]] != matrix(centsep[[dist_type]][i,],
                                           nrow=nrow(x),
                                           ncol=dist_cols,
                                           byrow=T)
          ifelse(is.na(d), 1, d)
        }, simplify = 'array')
        if(k == 1) {
          hlp[[dist_type]] <- array(hlp[[dist_type]],
                                    dim=c(nrow(x), dist_cols, k),
                                    dimnames=list(NULL,
                                                  names(which(dist_type==genDist)), NULL))
        }
      }
      
      if(dist_type == 'distJaccard') {
        
        centrep <- sapply(1:k, \(i) {
          matrix(centsep[[dist_type]][i,],
                 nrow=nrow(x),
                 ncol=dist_cols,
                 byrow=T)
        }, simplify = 'array')
        
        hlp[[dist_type]] <- sapply(1:k, \(i) {
          d <- xsep[[dist_type]] != centrep[,,i]
          d <- ifelse(is.na(d), 1, d) #das sollte im nächsten Schritt eh auch
          #mitgehandelt werden, aber zur Sicherheit explizit
          d[xsep[[dist_type]]==0] <- 1 #außerdem könnte ich mir den Schritt eh
          #sparen, die Info muss ja sowieso ins delta
          d
        }, simplify = 'array')
        if(k == 1) {
          hlp[[dist_type]] <- array(hlp[[dist_type]],
                                    dim=c(nrow(x), dist_cols, k),
                                    dimnames=list(NULL,
                                                  names(which(dist_type==genDist)), NULL))
        }
        
        delta <- sapply(1:k, \(i) {
          dlt <- xsep[[dist_type]]+centrep[,,i]
          dlt <- ifelse(is.na(dlt), 0, dlt)
          tmp_dlt <- array(delta[,,i,drop=F],
                           dim=dim(x), dimnames=dimnames(x))
          tmp_dlt[,genDist==dist_type][dlt==0] <- 0
          tmp_dlt
        }, simplify = 'array')
      }
    }
    
    hlp <- sapply(1:k, \(i) {
      # arr <- lapply(hlp, \(y) y[,,i, drop=F])
      # dims <- sapply(lapply(arr, dim), `[`, 2) #equivalent to sapply(hlp, ncol)
      # 
      # arr[dims==1] <- lapply(arr[dims==1],
      #                        \(y) matrix(y, dimnames=dimnames(y)[1:2]))
      # arr[dims>1] <- lapply(arr[dims>1], \(y) y[,,1])
      
      arr <- sapply(hlp, \(y) {
        arr <- y[,,i, drop=F]
        array(arr, dim=dim(arr)[1:2], dimnames=dimnames(arr)[1:2])
      })
      
      do.call(cbind, arr)
    }, simplify = 'array')
    hlp <- hlp[,colnames(x),,drop=F] #ordering by original column order
    
    z <- sapply(1:k,
                \(i) rowSums(hlp[,,i]*delta[,,i])/rowSums(delta[,,i]))
    return(z)
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
  #Error in y[, , i, drop = F] : incorrect number of dimensions
  kcca(dat, k, family=kccaFamilyGower(xmethods=.ChooseVarDists(sapply(dat, data.class))))
  #same error
  kcca(datNA, k, family=kccaFamilyGower())
  #same error
  
  #Error occurs within the allcent step
  #so:
  TestFam <- kccaFamilyGower()
  allcentenvir <- environment(TestFam@allcent)
  TestFam@allcent
  TestFam@allcent <- function(x, cluster, k=max(cluster, na.rm=TRUE))
  {
    centers <- matrix(NA, nrow=k, ncol=ncol(x))
    for(n in 1:k){
      browser()
      if(sum(cluster==n, na.rm=TRUE)>0){
        centers[n,] <- z@cent(x[cluster==n,,drop=FALSE])
      }
    }
    centers
  }
  environment(TestFam@allcent) <- allcentenvir
  kcca(dat,k,TestFam)
  #error occurs only for n==2
  ##cluster 2 occurs only once, all the others occurr multiple times
  # Browse[1]> cluster
  # [1] 2 3 3 1 1 3 3 3 4 4
  # Browse[1]> x[cluster==1,,drop=F]
  # cont   bin_sym  bin_asym ord_levmis ord_levfull       nom
  # [1,] 0.1860465 0.1162791 0.1162791  0.5813953   0.1162791 0.4651163
  # [2,] 0.7209302 0.1162791 0.1162791  0.2325581   0.4651163 0.3488372
  # Browse[1]> x[cluster==2,,drop=F]
  # cont bin_sym  bin_asym ord_levmis ord_levfull       nom
  # [1,]    1       0 0.1162791  0.2325581   0.2325581 0.1162791
  #soooo wait, would the kcca run if I repeated dat?
  kcca(dat[rep(1:nrow(dat), 5),], k, kccaFamilyGower())
  #YAAAS! so the issue doesn't 'really' lie in my code, but rather cause
  #the example is so small
  #well, gonna fix this anyway, but this is good news
  allcentenvir$z@cent
  allcentenvir$z@cent(x[1:3,,drop=F])
  allcentenvir$z@cent(x[1,,drop=F])
  #let's see whether the error occurs in the normal 'cent' as well,
  #cause this eval(bquote...) thing is a bit hard to debug
  #need to add genDist to the environment of @cent, cause this one
  #wasn't 'primed' in kcca:
  environment(TestFam@cent)$genDist <- .ChooseVarDists(dat)
  TestFam@cent(x[1:3,,drop=F])
  TestFam@cent(x[1,,drop=F]) #grreeeat, error created. Now let's debug
  centOptimNA <- function(x, dist) function(x, dist) {
    browser()
    foo <- function(p)
      sum(dist(x, matrix(p, nrow=1)), na.rm=TRUE)
    optim(colMeans(x, na.rm=TRUE), foo)$par
  }
  TestFam@cent(x[1,,drop=F]) #noup, doesn't access it
  rm(centOptimNA)
  h <- datmatrix[1,,drop=F]
  foo <- function(p) sum(distGower(h, matrix(p, nrow=1),
                                   genDist=.ChooseVarDists(dat)), na.rm=T)
  optim(colMeans(h, na.rm=T), foo)$par #Error occurs
  foo(colMeans(h, na.rm=T)) #Error occurs
  distGower(h, matrix(colMeans(h, na.rm=T), nrow=1),
            genDist=.ChooseVarDists(dat)) |> sum(na.rm=T) #Error occurs
  distGower(h, matrix(colMeans(h, na.rm=T), nrow=1),
            genDist=.ChooseVarDists(dat)) #Error occurs
  #so, let's debug within distGower
  rm(foo)
  #in distGower, the line y[,,i,drop=F] is found only in line88
  #hi
  #on Monday: run distGower with h; and also with datmatrix instead of h,
  #so I know what it's actually supposed to look like
  #go on debugging there
  
  
  #hi: right now: seems to work on distGower(datmatrix, centnoNAscaled, .ChooseVarDists(dat));
  #and distGower(x, centers, .ChooseVarDists(x)).
  #test these some more, and try out other edge cases
  #specifically with nrow(x)+nrow(centers) == 1
}
