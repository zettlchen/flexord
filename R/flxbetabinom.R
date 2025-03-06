# Based on Ivan Kondofersky's code from his bachelor's thesis.
# Added some light regularization


lbeta1 <- function(x, size, a, b) {
    s <- seq(from=0, to=size, by=1)
    uniquelb <- lbeta(a+s, b+size-s)
    res1 <- uniquelb[x+1L]
    res1
}

digamma1 <- function(x, a, size) {
    uniquedg <- digamma(seq(from=0, to=size, by=1) + a)
    res1 <- uniquedg[x+1L]
    res1
}


dbetabinom <- function(x, size, a, b, log=FALSE) {
    z <- lbeta1(x, size, a, b) - lbeta(a,b)
    if(!any(is.finite(z))) {
        warning("z not finite")
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

    c(sum((w+a2)*grad1), sum((w+a2)*grad2))
}

BBmle <- function(x, size=NULL, w=1, alpha2=0, eps=sqrt(.Machine$double.eps)) {
    N <- ncol(x)
    if (is.null(size)) {
        size <- apply(x,2,max, na.rm=TRUE)
    } else {
        size <- rep(size,length.out=N)
    }
    res <- matrix(NA_real_, nrow=2, ncol=N)
    for(i in seq_len(N)){
        res[,i] <- optim(c(1,1), fn=BBlogLikReg, gr=BBlogLikGradReg,
                         x=x[,i], size=size[i], w=w,
                         alpha2=alpha2,
                         control=list(fnscale=-1),
                         method="L-BFGS-B", lower=c(eps, eps))$par
    }
    rownames(res) <- c("alpha", "beta")
    colnames(res) <- colnames(x)
    res
}

#' FlexMix Driver for Regularized Beta-Binomial Mixtures
#'
#' This model driver can be used to cluster data using the beta-binomial
#' distribution.
#' 
#' Using a regularization parameter `alpha2` greater than zero can be
#' viewed as adding `alpha2` observations equal to the population mean
#' to each component. This can be used to avoid degenerate solutions
#' (i.e., probabilites of 0 or 1). It also has the effect that
#' clusters become more similar to each other the larger `alpha2` is
#' chosen. For small values this effect is, however, mostly
#' negligible.
#'
#' @param size Number of trials (one or more).
#' @param alpha2 A non-negative scalar acting as regularization
#'     parameter. Can be regarded as adding `alpha2` observations
#'     equal to the population mean to each component.
#' @param eps Lower threshold for the shape parameters a and b.
#' @param formula A formula which is interpreted relative to the
#'     formula specified in the call to [flexmix::flexmix()] using
#'     [stats::update.formula()]. Only the left-hand side (response)
#'     of the formula is used. Default is to use the original model
#'     formula specified in [flexmix::flexmix()].
#' @return an object of class `"FLXC"`
#' @export
#' @references
#'
#' Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).  *Ordinal
#' Clustering with the flex-Scheme.* Austrian Journal of
#' Statistics. _Submitted manuscript_.
#'
#' Kondofersky, I (2008). *Modellbasiertes Clustern mit der
#' Beta-Binomialverteilung.* Bachelor's thesis,
#' Ludwig-Maximilians-Universität München.
#' 
#' @export
#' @example examples/betabinom.R
FLXMCbetabinom = function(formula=.~., size, alpha2=0, eps=sqrt(.Machine$double.eps)) {
    z <- new("FLXMC", weighted=TRUE, formula=formula, dist="mvbetabinom",
             name="model based beta-binomial clustering")

    stopifnot(is.numeric(eps), length(eps) == 1, eps >= 0)
    stopifnot(is.numeric(alpha2), length(alpha2) == 1, alpha2 >= 0)
    size <- as.integer(size)
    stopifnot(all(size >= 1))

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
        para <- list(size = rep(size, length=ncol(y)))
        para$ab = BBmle(y, size=para$size, w=w, alpha2=alpha2, eps=eps)

        para$prob <- apply(para$ab, 2, function(z) z[1]/sum(z))
        para$center <- para$prob * para$size
        with(para, eval(z@defineComponent))
    }

    z
}
