#' FlexMix Driver for Regularized Multivariate Normal Mixtures
#'
#' This model driver implements the regularization method as introduced by
#' Fraley and Raftery (2007) for multivariate normal mixtures. Default
#' parameters for the regularization are taken from that paper.
#' However we only implement the special case where the covariance matrix
#' is diagonal and different variance per variable. For more general applications
#' of normal mixtures see package \pkg{Mclust}.
#'
#' For the regularization the conjugate prior distributions for the normal
#' distirbution are used, which are:
#' * Normal prior with parameter `mu_p` and `sigma^2/kappa_p` for the mean
#' * Inverse Gamma prior with parameters `nu_p/2` and `xi_p^2/2` tor the
#'   variance
#'
#'  `mu_p` is computed from the data as the overall means across all components.
#'
#' A value for the scale hyperparameter `xi_p` may be specified directly.
#' Otherwise the empirical variance divided by the square of the number of
#' components is used as per Fraley and Raftery (2007). In which case the
#' number of components (parameter `G`) needs to be specified.
#' 
#'
#' @param formula A formula describing the normal components
#' @param G Number of components in the mixture model (not used if xi_p is given)
#' @param kappa_p Shrinkage parameter. Functions as if you added
#'                `kappa_p` observations according to the population mean to
#'                each component (hyperparameter for IG prior)
#' @param nu_p Degress of freedom (hyperparameter for IG prior)
#' @param xi_p Scale (hyperparameter for IG prior). If not given the empirical
#'             variance divided by the square of the number of components
#'             is used as per Fraley and Raftery (2007).
#' @importFrom methods new
#' @importFrom stats cov.wt
#' @importFrom mvtnorm dmvnorm
#' @import flexmix
#' @export
#' @return an object of class `"FLXC"`
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Fraley, C, Raftery AE (2007)
#'   *Bayesian Regularization for Normal Mixture Estimation and Model-Based Clustering.*
#'   Journal of Classification, 24(2), 155-181
#' @example examples/regnorm.R
FLXMCregnorm <- function(formula=.~., xi_p=NULL, kappa_p=0.01, nu_p=3, G=NULL) {
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregnorm")


    if(is.null(xi_p) && is.null(G)) {
        stop("either parameter xi_p or G is needed")
    }

    z@defineComponent <- function(para) {
        predict <- function(x, ...){
            matrix(para$center, nrow=nrow(x), ncol=length(para$center), byrow=TRUE)
        }

        logLik <- function(x, y) {
            mvtnorm::dmvnorm(y, mean=para$center, sigma=diag(para$s2), log=TRUE)
        }

        new("FLXcomponent",
            parameters=list(center=para$center, s2=para$s2,
                            mu_p = para$mu_p, var_data = para$var_data),
            logLik=logLik, df=para$df, 
            predict=predict)
    }

    z@fit <- function(x, y, w, component=NULL) {
        # From Fraley/Raftery
        n = nrow(y)

        if(length(component) == 0L) {
            component$mu_p = colMeans(y)
            #component$var_data = diag(cov.wt(y, method="unbiased")$cov)
            component$var_data = 1/(n-1) * vapply(seq_len(ncol(y)), \(i) {
                sum((y[,i] - component$mu_p[i])^2)
            }, double(1))
        }

        if(is.null(xi_p) && !is.null(G)) {
            xi_p2 = component$var_data / G^2
        } else {
            xi_p2 = xi_p
        }


        nk = sum(w)
        ykbar = colSums(w*y)/nk

        muhat1 = (nk*ykbar + kappa_p*component$mu_p)/(kappa_p + nk)

        s2hat_numer1 = xi_p2 + (kappa_p*nk)/(kappa_p+nk)*(ykbar - component$mu_p)^2
        s2hat_numer2 = rowSums(w * (t(y) - ykbar)^2)
        s2hat_denom = nu_p + nk + 3
        s2hat = (s2hat_numer1+s2hat_numer2) / s2hat_denom


        para = list(center = muhat1, s2 = s2hat, df = 2*ncol(y))
        z@defineComponent(c(para, component))
    }

    z
}


