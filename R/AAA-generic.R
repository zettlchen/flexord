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
  if (is.null(xrange) || identical(xrange, 'all')) {
    rng <- function(x) {
      rep(range(x, na.rm=TRUE), ncol(x)) |>
          matrix(nrow=2)
    }
  } else if (identical(xrange, 'columnwise')) {
    rng <- function(x) {
        apply(x, 2, range, na.rm=TRUE)
    }
  } else if (is.vector(xrange, mode='numeric')) {
    rng <- function(x) {
        rep(xrange, ncol(x)) |>
            matrix(nrow=2)
    }
  } else if (is.list(xrange)) {
    rng <- function(x) {
        if (length(xrange) != ncol(x))
            stop('Either supply 1 range vector, or list of ranges for all variables')
        unlist(xrange) |> matrix(nrow=2)
    } 
  } else {
      stop("xrange not correctly specified")
  }
  return(rng)
}

