#Stand 17.01.25

#storage file for helpers that are used in more than 1 function,
#and archived functions

#helper to handle the different options in the xrange parameter.
#just prepares the function that will be used on x, the resulting
#function is applied within kcca
#used in: centMin (centroidFunctions.R); kccaFamilyGower (distGower.R)
# @param xrange: response level range of the variables. Implemented options:
#        'all': range of all x, over all variables
#        'columnwise': columnwise range of x
#        c(lower, upper): range vector specified by user, upper and lower
#                         limit for all of x
#        list(x1= c(lower, upper), x2=c(lower, upper), ...): list with user
#                 specified range vectors for each variable to be scaled.
.rangeMatrix <- function(xrange) {
  if(all(xrange=='all')) {
    rng <- function(x) {
      rep(range(x, na.rm=T), ncol(x)) |>
        matrix(nrow=2)
    }
  } else if(all(xrange=='columnwise')) {
    rng <- function(x) {
      apply(x, 2, range, na.rm=T)
    }
  } else if(is.vector(xrange, mode='numeric')) {
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

#archived scaling option after Gower that scales over all variable types.
#Easier to code (cause I don't need to extract xclass from the parent frame),
#and works in cases where xrange is variable specific ('columnwise', or list of range vectors).
#However, in the case of data wide ranges, it destroys the binary variables.
# @param x a numerically coded matrix.
# @param rangeMatrix expects a function that has been previously created
#                    with `.rangeMatrix(xrange)`. If it is NULL, it is created on x.
# @param xrange is a compatibility parameter so the helper runs outside
#               of the `kccaFamilyGower` concept, but within `kccaFamilyGower`,
#               `.rangeMatrix(xrange)` is run previously
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