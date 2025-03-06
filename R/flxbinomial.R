#' FlexMix Driver for Regularized Binomial Mixtures
#'
#' This model driver can be used to cluster data using the binomial
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
#' Parameter estimation is achieved using the MAP estimator for each
#' component and variable using a Beta prior.
#' 
#' @param size Number of trials (one or more).
#' @param alpha2 Regularization parameter. Can be regarded as adding
#'     `alpha2` observations equal to the population mean to each
#'     component.
#' @param eps A numeric value in [0, 1). When greater than zero,
#'     probabilities are truncated to be within in \[eps, 1-eps\].
#' @param formula A formula which is interpreted relative to the
#'     formula specified in the call to [flexmix::flexmix()] using
#'     [stats::update.formula()]. Only the left-hand side (response)
#'     of the formula is used. Default is to use the original model
#'     formula specified in [flexmix::flexmix()].
#' @return an object of class `"FLXC"`
#' @export
#' @references - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B
#'     (2025).  *Ordinal Clustering with the flex-Scheme.* Austrian
#'     Journal of Statistics. _Submitted manuscript_.
#' @example examples/binomial.R
FLXMCbinomial = function(formula=.~., size = NULL, alpha2=0, eps=0)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCbinomial")

    stopifnot(is.numeric(eps), length(eps) == 1, eps >= 0, eps < 1)
    stopifnot(is.numeric(alpha2), all(alpha2 >= 0))
    
    z@preproc.y <- function(y) {
        if(any(y < 0, na.rm=TRUE))
            stop("negative values are not allowed for the binomial family")
        y
    }

    defineComponent <- function(component, probs, df) {
        predict <- function(x, ...) {
            stop("not implemented")
        }

        logLik <- function(x, y) {
            probs = component$probs
            ty <- t(y)
            bc <- lchoose(size, ty)

            l <- bc + ty*log(probs) + (size-ty)*log(1-probs)
            cs <- colSums(l, na.rm=TRUE)
            if(any(!is.finite(cs))) stop("non-finite values while calculating log-likelihood")
            cs
        }

        new("FLXcomponent",
            parameters=component,
            logLik=logLik, predict=predict,
            df=component$df)
    }

    z@fit <- function(x, y, w, component) {
        if(length(component) == 0) {
            component$has_na <- anyNA(y)
            component$which_na <- which(is.na(y))

            if(is.null(size)) {
                component$size <- apply(y, 2, max, na.rm=TRUE)
            } else {
                component$size = size
            }
            component$ymarg <- colMeans(y, na.rm=TRUE)/component$size
            component$b_alpha <- component$ymarg*alpha2
            component$b_beta <- (1-component$ymarg)*alpha2
        }

        if(component$has_na) {
            p <- with(component,
                      (b_alpha + colSums(w*y, na.rm=TRUE)) /
                      (b_alpha+b_beta+size*colSums(w, na.rm=TRUE)))
        } else {
            p <- with(component,
                      (b_alpha + colSums(w*y, na.rm=TRUE)) /
                      (b_alpha+b_beta+size*sum(w, na.rm=TRUE)))
        }

        if(eps > 0) {
            p <- ifelse(p >= 1-eps, 1-eps, ifelse(p <= eps, eps, p))
        }

        component$probs = p
        component$df = ncol(y)
        defineComponent(component)
    }

    z
}


