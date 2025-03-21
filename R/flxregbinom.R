#' FlexMix Driver for Regularized Binomial Mixtures
#'
#' This model driver can be used to cluster data using the binomial
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
#' Parameter estimation is achieved using the MAP estimator for each
#' component and variable using a Beta prior.
#' 
#' @param size Number of trials (one or more).
#' @param alpha A non-negative scalar acting as regularization
#'     parameter. Can be regarded as adding
#'     `alpha` observations equal to the population mean to each
#'     component.
#' @param eps A numeric value in [0, 1). When greater than zero,
#'     probabilities are truncated to be within in \[eps, 1-eps\].
#' @param formula A formula which is interpreted relative to the
#'     formula specified in the call to [flexmix::flexmix()] using
#'     [stats::update.formula()]. Only the left-hand side (response)
#'     of the formula is used. Default is to use the original model
#'     formula specified in [flexmix::flexmix()].
#' @param has_na Boolean whether the data set may contain NA values. Default
#'               is FALSE. For data sets without NAs this parameter does not
#'               influence the estimates except that it's slightly faster
#'               when the absence of NAs can be assumed.
#' @return an object of class `"FLXC"`
#' @export
#' @references - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B
#'     (2025).  *Ordinal Clustering with the flex-Scheme.* Austrian
#'     Journal of Statistics. _Submitted manuscript_.
#' @example examples/binomial.R
#' @importFrom stats dbinom
FLXMCregbinom = function(formula=.~., size = NULL, has_na=FALSE, alpha=0, eps=0)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregbinom")

    stopifnot(is.numeric(eps), length(eps) == 1, eps >= 0, eps < 1)
    stopifnot(is.numeric(alpha), length(alpha) == 1, alpha >= 0)
    
    z@preproc.y <- function(y) {
        if(any(y < 0, na.rm=TRUE))
            stop("negative values are not allowed for the binomial family")
        y
    }

    defineComponent <- function(component) {
        predict <- function(x, ...) {
            matrix(component$probs * component$size, nrow = nrow(x), ncol = length(component$probs),
                   byrow = TRUE)
        }

        logLik <- function(x, y) {
            llh <- dbinom(t(y), size = component$size, prob = component$probs, log = TRUE)
            colSums(llh, na.rm = TRUE)
        }

        new("FLXcomponent",
            parameters=component,
            logLik=logLik, predict=predict,
            df=component$df)
    }

    z@fit <- function(x, y, w, component) {
        if(length(component) == 0) {
            #component$has_na <- anyNA(y)

            if(is.null(size)) {
                component$size <- apply(y, 2, max, na.rm=TRUE)
            } else {
                component$size = size
            }
            component$ymarg <- colMeans(y, na.rm=TRUE)/component$size
            component$b_alpha <- component$ymarg*alpha
            component$b_beta <- (1-component$ymarg)*alpha
        }

        if(has_na) {
            p <- with(component,
                      (b_alpha + colSums(w*y, na.rm=TRUE)) /
                      (b_alpha+b_beta+size*colSums(w * !is.na(y))))
        } else {
            p <- with(component,
                      (b_alpha + colSums(w*y)) /
                      (b_alpha+b_beta+size*sum(w)))
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


