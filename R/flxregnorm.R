#' FlexMix driver for regularized multivariate normal mixtures
#'
#' This model implements the regularization method as introduced in
#' Fraley & Raftery (...) for multivariate normal mixtures.
#' The covariance matrix for each component is assumed to be diagonal.
#' @param formula A formula describing the normal components
#' @param G Number of components in the mixture model
#' @param kappa_p Regularization parameter. Functions as if you added
#'                kappa_p observations according to the population mean to the
#'                data
#' @param nu_p Regularization parameter.
#' 
#' @importFrom methods new
#' @importFrom stats cov.wt
#' @importFrom mvtnorm dmvnorm
#' 
#' @export
FLXMCregnorm <- function(formula=.~., G, kappa_p=0.01, nu_p=3)
{
    z <- methods::new("FLXMC", weighted=TRUE, formula=formula,
                      name="FLXMCregnorm")

    force(G)

    z@defineComponent <- function(para) {
        predict <- function(x, ...){
            matrix(para$center, nrow=nrow(x), ncol=length(para$center), byrow=TRUE)
        }

        logLik <- function(x, y) {
            mvtnorm::dmvnorm(y, mean=para$center, sigma=diag(para$s2), log=TRUE)
        }

        methods::new("FLXcomponent",
                     parameters=list(center=para$center, s2=para$s2),
                     logLik=logLik, df=para$df, 
                     predict=predict)
    }

    z@fit <- function(x, y, w, ...) {
        # From Fraley/Raftery
        n = nrow(y)
        mu_p = colMeans(y)

        var_data = stats::cov.wt(y, method="unbiased")$cov |> diag()
        xi_p2 = var_data / G^2


        nk = sum(w)
        ykbar = colSums(w*y)/nk

        muhat1 = (nk*ykbar + kappa_p*mu_p)/(kappa_p + nk)


        s2hat_numer1 = xi_p2 + (kappa_p*nk)/(kappa_p+nk)*(ykbar - mu_p)^2
        s2hat_numer2 = rowSums(w * (t(y) - ykbar)^2)
        s2hat_denom = nu_p + nk + 3
        s2hat = (s2hat_numer1+s2hat_numer2) / s2hat_denom


        para = list(center = muhat1, s2 = s2hat, df = 2*ncol(y))
        z@defineComponent(para)
    }

    z
}


