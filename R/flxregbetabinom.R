# Based on Ivan Kondofersky's code from his bachelor's thesis.
# Added some light regularization


lbeta1 <- function(x, size, a, b) {
    s <- 0:size
    uniquelb <- lbeta(a+s, b+size-s)
    uniquelb[x+1L]
}

digamma1 <- function(x, a, size) {
    uniquedg <- digamma((0:size) + a)
    uniquedg[x+1L]
}


dbetabinom <- function(x, size, a, b, log=FALSE) {
    z <- lbeta1(x, size, a, b) - lbeta(a,b)
    if (log) z else exp(z)
}

BBlogLikGrad <- function(ab, x, size, w=1) {
    c(sum(w*(digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))),
             
      sum(w*(digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))))
}

BBlogLikReg <- function(ab, x, size, w=1, alpha=0) {
    a2 = alpha/length(x)
    dens = dbetabinom(x, size, ab[1], ab[2], log=TRUE)
    sum((w+a2)*dens)
}

BBlogLikGradReg <- function(ab, x, size, w=1, alpha=0) {
    grad1 = (digamma1(x, ab[1], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[1]) + digamma(ab[1]+ab[2]))
    grad2 = (digamma1(size-x, ab[2], size=size) -
             digamma(size+ab[1]+ab[2]) -
             digamma(ab[2]) + digamma(ab[1]+ab[2]))

    a2 = alpha/length(x)

    c(sum((w+a2)*grad1), sum((w+a2)*grad2))
}

BBmle <- function(x, size=NULL, w=1, alpha=0, eps=sqrt(.Machine$double.eps)) {
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
                         alpha=alpha,
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
#' Using a regularization parameter `alpha` greater than zero can be
#' viewed as adding `alpha` observations equal to the population mean
#' to each component. This can be used to avoid degenerate solutions
#' (i.e., probabilites of 0 or 1). It also has the effect that
#' clusters become more similar to each other the larger `alpha` is
#' chosen. For small values this effect is, however, mostly
#' negligible.
#'
#' @param size Number of trials (one or more). Default `NULL` implies
#'     that the number of trials is inferred columnwise by the
#'     maximum value observed.
#' @param alpha A non-negative scalar acting as regularization
#'     parameter. Can be regarded as adding `alpha` observations
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
FLXMCregbetabinom = function(formula=.~., size=NULL, alpha=0, eps=sqrt(.Machine$double.eps)) {
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregbetabinom")

    stopifnot(is.numeric(eps), length(eps) == 1, eps >= 0)
    stopifnot(is.numeric(alpha), length(alpha) == 1, alpha >= 0)

    z@preproc.y <- function(y) {
        if (anyNA(y)) 
            stop("NAs are not allowed")
        if (any(y < 0))
            stop("negative values are not allowed for the binomial family")
        if (is.null(size)) {
            size <- apply(y, 2, max)
        } else {
            size <- as.integer(size)
            size <- rep(size, length = ncol(y))
            if (any(apply(y, 2, max) > size)) 
                stop("values larger than size not allowed (values need to be in 0:size)")
        }
        stopifnot(size >= 1)
        
        attr(y, "size") <- size
        y
    }

    z@defineComponent <- function(para) {
        predict <- function(x, ...) {
            matrix(para$prob * para$size,
                   nrow=nrow(x), ncol=length(para$prob), byrow=TRUE)
        }

        logLik <- function(x, y) {
            z = lapply(seq_len(ncol(y)), \(j) {
                dbetabinom(y[,j], size=para$size[j], a=para$ab[1,j], b=para$ab[2,j], log=TRUE)
            }) |> do.call(cbind, args=_)
            rowSums(z, na.rm=TRUE)
        }
        
        new("FLXcomponent", 
            parameters = list(a = para$ab[1, ], b = para$ab[2, ]),
            logLik = logLik, predict = predict, df = length(para$ab))
    }

    z@fit <- function(x, y, w) {
        para <- list(size = attr(y, "size"))
        para$ab <- BBmle(y, size=para$size, w=w, alpha=alpha, eps=eps)
        para$prob <- apply(para$ab, 2, function(z) z[1]/sum(z))

        z@defineComponent(para = para)
    }

    z
}
