# 06.12.24

#script for GDM2 distance and its helper, main documentation block in distSimMatch.R

#' @rdname distances
#' @importFrom stats ecdf
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

# @param x a numerically coded matrix.
# @param rangeMatrix expects a function that has been previously created
#                    with `.rangeMatrix(xrange)`. If it is NULL, it is created on x.
# @param xrange is a compatibility parameter so the helper runs outside
#               of the `kccaFamilyGDM2` concept, but within `kccaFamilyGDM2`,
#               `.rangeMatrix(xrange)` is run previously
.projectIntofx <- function(x, rangeMatrix=NULL,
                           xrange=NULL){
  
  if(is.null(rangeMatrix)) {
    rng <- .rangeMatrix(xrange)
  } else {
    rng <- rangeMatrix
  }
  
  rng <- rng(x)
  
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
    list(epdf=epdf, ecdf=stats::ecdf(x[,y]))
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

