#' FlexMix Driver for Regularized Multivariate Normal Mixtures
#'
#' This model driver implements the regularization method as
#' introduced by Fraley and Raftery (2007) for univariate normal
#' mixtures. Default parameters for the regularization according to
#' that paper may be obtained using `FLXMCregnorm_defaults()`. We
#' extend this to the multivariate case assuming independence between
#' variables within components, i.e., we only implement the special
#' case where the covariance matrix is diagonal. For more general
#' applications of normal mixtures see package \pkg{mclust}.
#'
#' For the regularization the conjugate prior distributions for the
#' normal distribution are used, which are:
#' 
#' * Normal prior with parameter `mu_p` and `sigma^2/kappa_p` for the mean.
#' * Inverse Gamma prior with parameters `nu_p/2` and `zeta_p^2/2` for the
#'   variance.
#'
#' @param formula A formula which is interpreted relative to the formula
#'        specified in the call to [flexmix::flexmix()] using
#'        [stats::update.formula()]. Only the left-hand side (response)
#'        of the formula is used. Default is to
#'        use the original model formula specified in [flexmix::flexmix()].
#' @param params Prior parameters for normal mixtures. You may obtain 
#'        default values according to Fraley and Raftery (2007) using 
#'        `FLXMCregnorm_defaults()`.
#'        As the prior depends on the number of components it is probably not
#'        advisable to run `stepFlexmix` with more than one value of `k` at
#'        a time.
#' @importFrom methods new
#' @importFrom mvtnorm dmvnorm
#' @import flexmix
#' @export
#' @return An object of class `"FLXC"`.
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).
#'   *Ordinal Clustering with the flex-Scheme.*
#'   Austrian Journal of Statistics. _Submitted manuscript_.
#' - Fraley, C, Raftery, AE (2007)
#'   *Bayesian Regularization for Normal Mixture Estimation and Model-Based Clustering.*
#'   Journal of Classification, 24(2), 155-181
#' @seealso FLXMCregnorm_defaults
#' @example examples/regnorm.R
FLXMCregnorm <- function(formula=.~., params) {
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregnorm")

    z@defineComponent <- function(para) {
        predict <- function(x, ...){
            matrix(para$center, nrow=nrow(x), ncol=length(para$center), byrow=TRUE)
        }

        logLik <- function(x, y) {
            mvtnorm::dmvnorm(y, mean=para$center, sigma=diag(para$s2), log=TRUE)
        }

        new("FLXcomponent",
            parameters=list(center = para$center, s2 = para$s2),
            logLik=logLik, df=para$df, 
            predict=predict)
    }

    z@fit <- function(x, y, w, component=NULL) {
        nk <- sum(w)
        ykbar <- colSums(w*y)/nk

        muhat1 <- (nk*ykbar + params$kappa_p*params$mu_p)/(params$kappa_p + nk)

        s2hat_numer1 <- params$zeta_p +
            (params$kappa_p*nk)/(params$kappa_p+nk)*(ykbar - params$mu_p)^2
        s2hat_numer2 <- rowSums(w * (t(y) - ykbar)^2)
        s2hat_denom <- params$nu_p + nk + 3
        s2hat <- (s2hat_numer1 + s2hat_numer2) / s2hat_denom


        para = list(center = muhat1, s2 = s2hat, df = 2*ncol(y))
        z@defineComponent(c(para, component))
    }

    z
}


#' Data-Driven Default Parameters for Regularized Normal Mixtures 
#'
#' Determines the default values for regularized univariate normal
#' mixtures as proposed by Fraley and Raftery (2007) based on the data
#' set to be clustered and the number of components in the mixture
#' model.
#' 
#' @param x The data set to be clustered. Should be the same data set
#'     as is used in [flexmix::flexmix()]'s model formula.
#' @param k Number of components assumed for the mixture model (not
#'     used if `zeta_p` is given).
#' @param kappa_p Shrinkage parameter. Corresponds to adding `kappa_p`
#'     observations according to the population mean to each component
#'     (hyperparameter for IG prior)
#' @param nu_p Degress of freedom (hyperparameter for IG prior)
#' @param zeta_p Scale (hyperparameter for IG prior). If not given the
#'     empirical variance divided by the square of the number of
#'     components is used as per Fraley and Raftery (2007).
#'
#' `mu_p` is computed from the data as the overall means across all
#'     components.
#'
#' A value for the scale hyperparameter `zeta_p` may be specified directly.
#' Otherwise the empirical variance divided by the square of the number of
#' components is used as per Fraley and Raftery (2007). In which case the
#' number of components (parameter `k`) needs to be specified.
#' @return A named list with values for `mu_p`, `kappa_p`, `nu_p` and `zeta_p`.
#' 
#' @export
FLXMCregnorm_defaults <- function(x, zeta_p=NULL, kappa_p=0.01, nu_p=3, k=NULL) {
    if (is.null(zeta_p) && is.null(k)) {
        stop("either parameter zeta_p or k is needed")
    }
    stopifnot(is.numeric(kappa_p), length(kappa_p) == 1)
    stopifnot(is.numeric(nu_p), length(nu_p) == 1)

    params <- list(mu_p = colMeans(x),
                   kappa_p = kappa_p,
                   nu_p = nu_p)

    if (is.null(zeta_p)) {
        stopifnot(is.numeric(k), length(k) == 1)
        stopifnot(is.matrix(x) || is.data.frame(x))
        k <- as.integer(k)
        n <- nrow(x)
        var_data <- 1/(n-1) * vapply(seq_len(ncol(x)), \(i) {
            sum((x[,i] - params$mu_p[i])^2)
        }, double(1))
        zeta_p <- var_data / k^2
    } 

    params$zeta_p <- zeta_p
    params
}



