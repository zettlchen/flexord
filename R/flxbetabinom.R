# Based on Ivan Kondofersky's code from his bachelor's thesis.
# Added some light regularization


lbeta1 <- function(x, size, a, b) {
    #unique(cbind(x, size-x))
    
    s <- seq(from=0, to=size, by=1)
    uniquelb <- lbeta(a+s, b+size-s)

    res1 <- uniquelb[x+1L]
    #res2 <- lbeta(a+x, b+size-x)

    #if(!all(res1==res2)) browser() else cat("ok\n")

    #sort(unique(res2))
    #sort(unique(uniquelb))

    res1
}

digamma1 <- function(x, a, size) {
    uniquedg <- digamma(seq(from=0, to=size, by=1) + a)
    res1 <- uniquedg[x+1L]
    #res2 <- digamma(x+a)

    #if(!all(res1==res2)) browser()

    # unique(abs(res1-res2))
    #sort(unique(digamma(x+px)))
    res1
}


dbetabinom <- function(x, size, a, b, log=FALSE) {
    #z <- lchoose(size, x) + lbeta(a+x, b+size-x) - lbeta(a,b)
    #z <- lbeta(a+x, b+size-x) - lbeta(a,b)
    z <- lbeta1(x, size, a, b) - lbeta(a,b)
    if(!any(is.finite(z))) {
        cat("z not finite\n")
        return(NA)
    }
    if(log) z else exp(z)
}

BBlogLikGrad <- function(ab, x, size, w=1) {
    c(sum(w*(digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))),
             
      sum(w*(digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))))
}

BBlogLikReg <- function(ab, x, size, w=1, alpha2=0) {
    a2 = alpha2/length(x)
    # bad
    if(FALSE) {
        sum(w * dbetabinom(x, size, ab[1], ab[2], log=TRUE)) +
            sum(a2 * dbetabinom(x, size, ab[1], ab[2], log=TRUE))
    }

    # slightly less bad
    if(FALSE) {
        dens = dbetabinom(x, size, ab[1], ab[2], log=TRUE)
        sum(w * dens) + sum(a2 * dens)
    }

    # bit better
    dens = dbetabinom(x, size, ab[1], ab[2], log=TRUE)
    sum((w+a2)*dens)
}

BBlogLikGradReg <- function(ab, x, size, w=1, alpha2=0) {
    grad1 = (digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))
    grad2 = (digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))

    a2 = alpha2/length(x)

    if(FALSE) {
        c(sum(w*grad1) + sum(a2*grad1),
          sum(w*grad2) + sum(a2*grad2))
    }

    c(sum((w+a2)*grad1), sum((w+a2)*grad2))
}

BBmle <- function(x, size=NULL, w=1, alpha2=0, eps=sqrt(.Machine$double.eps))
{
    N<-ncol(x)
    if (is.null(size)) size <- apply(x,2,max, na.rm=TRUE)
    else size<-rep(size,length.out=N)
    res<-matrix(NA, nrow=2, ncol=N)
    for(i in seq_len(N)){
        res[,i]<-optim(c(1,1), fn=BBlogLikReg, gr=BBlogLikGradReg,
                       x=x[,i], size=size[i], w=w,
                       alpha2=alpha2,
                       control=list(fnscale=-1),
                       method="L-BFGS-B", lower=c(eps, eps))$par
    }
    rownames(res)<-c("alpha","beta")
    colnames(res)<-colnames(x)
    res
}

#' FlexMix driver for regularized beta-binomial mixtures
#'
#' This model driver can be used to cluster data using the beta-binomial
#' distribution.
#' 
#' Using a regularization parameter `alpha2` greater than zero
#' acts as adding `alpha2` observations conforming to the population
#' mean to each component. This can be used to avoid degenerate
#' solutions. It also has the effect
#' that clusters become more similar to each other the larger
#' `alpha2` is chosen. For small values it is mostly negligible however.
#'
#' @param size number of trials (zero or more)
#' @param alpha2 Regularization parameter. Can be regarded the same as
#'  adding `alpha2` observations conforming to the population mean to each
#'  component.
#' @param eps Lower threshold for the shape parameters a and b
#' @return an object of class FLXC
#' @export
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Kondofersky, Ivan (2008).
#'   *Modellbasiertes Clustern mit der Beta-Binomialverteilung.*
#'   Bachelor's thesis, Ludwig-Maximilians-Universität München
#' @export
#' @example examples/betabinom.R
FLXMCbetabinom = function(formula=.~., size, alpha2=0, eps=sqrt(.Machine$double.eps)) {

    z <- new("FLXMC", weighted=TRUE, formula=formula, dist="mvbetabinom",
             name="model based beta-binomial clustering")

    size <- as.integer(size)
    

    z@defineComponent <- expression({
        logLik <- function(x,y) {
            z <- y
            for(k in 1:ncol(y)) {
                z[,k] <- dbetabinom(y[,k], size=size[k], a=ab[1,k], b=ab[2,k], 
                                    log=TRUE)
            }
            rowSums(z, na.rm=TRUE)
        }
        predict <- function(x, ...) {
            matrix(center, nrow=nrow(x), ncol=length(center), byrow=TRUE)
        }

        new("FLXcomponent", 
            parameters=list(a=ab[1,], b=ab[2,], size=size, prob=prob, center=center),
            df=length(ab), logLik=logLik, predict=predict)
    })

    z@fit <- function(x,y,w) {
        para <- list(size=rep(size, length=ncol(y)))
        para$ab = BBmle(y, size=para$size, w=w, alpha2=alpha2, eps=eps)

        para$prob <- apply(para$ab, 2, function(z) z[1]/sum(z))
        para$center <- para$prob * para$size
        with(para, eval(z@defineComponent))
    }

    z
}
