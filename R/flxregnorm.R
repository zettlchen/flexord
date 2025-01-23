#' FlexMix driver for regularized multivariate normal mixtures
#'
#' This model implements the regularization method as introduced in
#' Fraley & Raftery (...) for multivariate normal mixtures.
#' The covariance matrix for each component is assumed to be diagonal.
#' @param formula A formula describing the normal components
#' @param G Number of components in the mixture model (not used if xi_p is given)
#' @param kappa_p Regularization parameter. Functions as if you added
#'                kappa_p observations according to the population mean to the
#'                data
#' @param nu_p Regularization parameter.
#' 
#' @importFrom methods new
#' @importFrom stats cov.wt
#' @importFrom mvtnorm dmvnorm
#' @import flexmix
#' @export
#' @return an object of class FLXC
#' @export
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' - Fraley, C, Raftery AE (2007)
#'   *Bayesian Regularization for Normal Mixture Estimation and Model-Based Clustering.*
#'   Journal of Classification, 24(2), 155-181
FLXMCregnorm <- function(formula=.~., xi_p=NULL, kappa_p=0.01, nu_p=3, G=NULL) {
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


