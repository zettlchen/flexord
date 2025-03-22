# Stand: 18.01.25

#script for Gower's distance and its helpers, main documentation block in distSimMatch.R

#old distGower_ordinal that was used in paper (written just for ordinal, not mixed)
#can be found in the AJS paper folder, is now replaced here as the other options work properly

# .ChooseVarDists: helper that maps default distances to 4 different variable types.
#                  A different option is to specifically provide a character vector
#                  of length ncol(x) to kccaFamilyGower(xmethods) that specifies the
#                  distance for each variable. In the latter case, this helper is
#                  circumvented.
# @param xclass Character vector of length=ncol(x) of classes of variables in x,
#               for example as obtained by sapply(data.frame(x), data.class).
#                   Default variable specific methods will be mapped to each variable class:
#                     - 'numeric' or 'integer': Euclidean distance ('distEuclidean'),
#                                               expects data to be scaled previously (f.i. by
#                                               .ScaleVarSpecific) 
#                     - 'logical': Jaccard distance ('distJaccard')
#                     - 'ordered': Manhattan distance ('distManhattan'),
#                                               expects data to be scaled previously (f.i. by
#                                               .ScaleVarSpecific) 
#                     - 'factor' (i.e. categorical): Simple Matching Distance ('distSimMatch')
#                  The treatment of the different variables can thus influenced in two ways:
#                  First, by specifing the xmethods parameter in kccaFamilyGower(),
#                  Secondly by their coding. F.i. for a symmetric logical variable,
#                  symmetric treatment can be achieved by coding them either as
#                  categorical or numeric instead of logical.
# written after: Gower (1971), Kaufman & Rousseeuw (1990).
.ChooseVarDists <- function(xclass) {
  
  if (!is.null(dim(xclass))) {#compatibility option, where .ChooseVarDists is used outside of kcca, and directly on x
    if (is.data.frame(xclass)) {
      xclass <- sapply(xclass, data.class)
    } else {
      xclass <- rep("numeric", ncol(xclass))
    }
  }

  xclass <- factor(xclass, c("numeric", "logical", "ordered", "factor"))
  if (anyNA(xclass)) {
      stop("no default method implement for one of the classes")
  }
  c("distEuclidean", "distJaccard", "distManhattan", "distSimMatch")[as.integer(xclass)] |>
      stats::setNames(names(xclass))
}


#.ScaleVarSpecific: Variable scaling after Gower and Kaufman+Rousseeuw (center around
#the minimum value, and divide by range). Only performed on numeric
#or ordinal variables. In the latter case, the function expects the
#variables to be coded numerically ranging from `1:max(level)` in steps of 1.
# @param x a numerically coded matrix. (or, within stepFlexclust: a dataframe without character columns)
# @param xclass character vector of either:
#               - variable classes as obtained by `data.class`; or:
#               - proposed variable-specific default distances, f.i.
#                 as obtained by `.ChooseVarDists(x)`;
#               with `length(xclass)==ncol(x)`
# @param rangeMatrix expects a function that has been previously created
#                    with `.rangeMatrix(xrange).` If it is NULL, it will be
#                    created on `xrange`.
# @param xrange is a compatibility parameter so the helper runs outside
#               of the `kccaFamilyGower` concept, but within `kccaFamilyGower`,
#               `.rangeMatrix(xrange)` is run previously
.ScaleVarSpecific <- function(x, xclass, rangeMatrix = NULL, xrange = NULL) {
  
  if (is.null(rangeMatrix)) {
    rng <- .rangeMatrix(xrange)
  } else {
    rng <- rangeMatrix
  }
  
  if (is.data.frame(x)) { #this happens only within stepFlexclust, within kcca, data.matrix(x) has already been run
    x <- data.matrix(x)
  }
  
  if (is.null(colnames(x))) {
      colnames(x) <- 1:ncol(x) #don't need this here but later in the dist #I think it's actually obsolete now
  }
  
  cols2scl <- xclass %in% c('numeric', 'ordered', #compatibility with xclass=sapply(dat, data.class)
                            'distEuclidean', 'distManhattan') #compatibility with xmethods
  
  if (sum(cols2scl) > 0) {
    
    rng <- rng(x[, cols2scl, drop=FALSE])
    scl <- apply(rng, 2, diff)
    scl <- ifelse(scl==0, 1, scl)
    
    x[, cols2scl] <- scale(x[, cols2scl], center=rng[1,], scale=scl)
  }

  return(x)
}

#' @rdname distances
#' @export
distGower <- function(x, centers, genDist) {
#.distGower_mixedType: helper function to calculate Gower's distance on _mixed variable types_, and/or
#in the presence of missing values
# @param x a numerically coded matrix. If `distEuclidean` or
#          `distManhattan` are to be used on some columns, these
#          columns need to be scaled previously, f.i. by
#          `.ScaleVarSpecific(x)`.
# @param centers a numerically coded matrix that is compatible with,
#                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#                `nrow(x)>=nrow(centers)`.
# @param distances character vector of variable specific distances,
#                  f.i. as derived from `.ChooseVarDists(x)`. Length
#                  needs to be equal to `ncol(x)`.
# @return A numeric array of dim `nrow(x)`X`ncol(x)`X`nrow(centers)`
#         that in each `ijk` contains the distance between `x[i,j]`
#         and `centers[k,j]`. Weighting with `delta`, and summing up
#         over `j` happens in the next step.
# @details NA handling: As proposed by Gower, and Kaufman+Rousseeuw. NAs
#          are filled with placeholder values, and handled in the weighting
#          step (Pairs of `c(x[i,j], centers[k,j])` with missing values receive
#          weight 0, and values that are present in `c(x[i,], centers[k,])` are
#          upweighted.
.distGower_mixedType <- function(x, centers, distances) {
  
  z <- array(0, dim=c(dim(x), nrow(centers)),
             dimnames=c(dimnames(x), list(NULL)))
  K <- nrow(centers)
  distE <- distances %in% 'distEuclidean'
  distM <- distances %in% 'distManhattan'
  distS <- distances %in% 'distSimMatch'
  distJ <- distances %in% 'distJaccard'
  if (sum(distE, distM, distS, distJ) != ncol(x)) {
      stop('Specified distance(s) not implemented in flexord::distGower')
  }
  
  if (any(distE)) {
    xE <- x[, distE, drop = FALSE]
    cE <- centers[, distE, drop = FALSE]
    for (k in 1:K) {
      z[, distE, k] <- t(sqrt((t(xE) - cE[k,])^2))
    }
  }
  if (any(distM)) {
    xM <- x[, distM, drop = FALSE]
    cM <- centers[, distM, drop = FALSE]
    for (k in 1:K) {
      z[, distM, k] <- t(abs(t(xM)-cM[k,]))
    }
  }
  if (any(distS)) {
    xS <- x[, distS, drop = FALSE]
    cS <- centers[, distS, drop = FALSE]
    for (k in 1:K) {
      z[, distS, k] <- t(t(xS) != cS[k,])
    }
  }
  if (any(distJ)) {
    xJ <- x[, distJ, drop = FALSE]
    cJ <- centers[, distJ, drop = FALSE]
    for (k in 1:K) {
      z[, distJ, k] <- t((t(xJ) + cJ[k,]) < 2)
    }
  }

  z[is.na(z)] <- 1 #just a placeholder, NA cases are removed by weights==0
  z
}

#.distGower_singleTypeNoNAs: helper function to calculate Gower's distance on _single variable types_,
#and no missing values.
# @param x a numerically coded matrix. If `distEuclidean` or
#          `distManhattan` are to be used on some columns, these
#          columns need to be scaled previously, f.i. by
#          `.ScaleVarSpecific(x)`.
# @param centers a numerically coded matrix that is compatible with,
#                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#                `nrow(x)>=nrow(centers)`.
# @param distances character vector of distance to be used on all variables.
#                  Length needs to be equal to `ncol(x)`, but needs to be the
#                  same for all columns. Can f.i. be derived from `.ChooseVarDists(x)`.
# @return A numeric matrix of dim `nrow(x)`X`nrow(centers)`
#         that in each `ik` contains the distance between `x[i,]`
#         and `centers[k,]`.
# @details
# This helper is added for the speed increase it provides when using
# Gower's distance on all numeric or all ordered variables without
# missing values.
# 
# The difference between this option, and using `distEuclidean` or 
# `distManhattan` directly is that 1) when using these distances within
# Gower's distance, they are scaled previously as described by Gower and
# Kaufman+Rousseeuw, and 2) the are then divided by the number of columns
# (as all weights will be 1). This results in a distance that ranges from 0-1.
# 
# Due to compatibility reasons, this function can also handle Jaccard distance,
# #or Simple Matching Distance over all variables. However, in this case it
# brings no advantages. Thus, it is recommended to use `flexclust::distJaccard`
# or `distSimMatch` directly.
# 
# The checks for this helper to run (`length(unique(distances))==1` and
# `!anyNA(x)` happen outside of it.
.distGower_singleTypeNoNAs <- function(x, centers, distances) {

  dists <- unique(distances)
  p <- ncol(x)
  
  if (dists == 'distEuclidean') {
    dstfnc <- flexclust::distEuclidean
  } else if (dists == 'distManhattan') {
    dstfnc <- flexclust::distManhattan
  } else if (dists == 'distSimMatch') {
    dstfnc <- distSimMatch
    p <- 1 
  } else if (dists == 'distJaccard') {
    dstfnc <- flexclust::distJaccard
    p <- 1 
  } else {
      stop('Specified distance not implemented in flexord::distGower')
  }

  dstfnc(x, centers)/p
}

#.delta: helper function to calculate weights for Gower's distance
# @param x a numerically coded matrix.
# @param centers a numerically coded matrix that is compatible with,
#                or a subset of, `x`. `ncol(x)==ncol(centers)`, and
#                `nrow(x)>=nrow(centers)`.
# @param distances character vector of variable specific distances,
#                  f.i. as derived from `.ChooseVarDists(x)`. Length
#                  needs to be equal to `ncol(x)`.
# @return A logical array of dim `nrow(x)`X`ncol(x)`X`nrow(centers)`.
#         `TRUE` in `ijk` if `x[i,j]` and `centers[k,j]` are both non-
#         missing, and, in the case of logical variables calculated with
#         `distJaccard`, not more than one variables is equal to 0; `FALSE`
#         otherwise. 
.delta <- function(x, centers, distances) {
  K <- nrow(centers)
  delta <- sapply(1:K,
                  \(k) t(!(is.na(t(x)) | is.na(centers[k,]))),
                  simplify='array')
  distJ <- distances == 'distJaccard'
  
  if (any(distJ)) {
    xJ <- x[, distJ, drop = FALSE]
    cJ <- centers[, distJ, drop = FALSE]
    for (k in 1:K) {
      delta[, distJ, k] <- t((t(xJ) + cJ[k,])>0)
    }
    delta[is.na(delta)] <- FALSE
  }
  
  delta
  
}

  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  
  if (length(unique(genDist))==1 && !anyNA(x)) {
    .distGower_singleTypeNoNAs(x, centers, distances = genDist)
  } else {
    delta <- .delta(x, centers, distances=genDist)
    z <- .distGower_mixedType(x, centers, distances = genDist)
    z <- sapply(1:nrow(centers), \(k) {
        rowSums(z[,,k,drop = FALSE] * delta[,,k,drop=FALSE]) / rowSums(delta[,,k,drop = FALSE])})
    z[is.nan(z)] <- 0 # see distJaccard
    
    z
  }
}

