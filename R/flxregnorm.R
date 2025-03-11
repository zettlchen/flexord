#' FlexMix Driver for Regularized Multivariate Normal Mixtures
#'
#' This model driver implements the regularization method as
#' introduced by Fraley and Raftery (2007) for univariate normal
#' mixtures. Default parameters for the regularization are taken from
#' that paper.  We extend this to the multivariate case assuming
#' independence between variables within components, i.e., we only
#' implement the special case where the covariance matrix is
#' diagonal. For more general applications of normal mixtures see
#' package \pkg{mclust}.
#'
#' For the regularization the conjugate prior distributions for the
#' normal distribution are used, which are:
#' 
#' * Normal prior with parameter `mu_p` and `sigma^2/kappa_p` for the mean.
#' * Inverse Gamma prior with parameters `nu_p/2` and `zeta_p^2/2` for the
#'   variance.
#'
#' `mu_p` is computed from the data as the overall means across all components.
#'
#' A value for the scale hyperparameter `zeta_p` may be specified directly.
#' Otherwise the empirical variance divided by the square of the number of
#' components is used as per Fraley and Raftery (2007). In which case the
#' number of components (parameter `G`) needs to be specified.
#' 
#'
#' @param formula A formula which is interpreted relative to the formula
#'        specified in the call to [flexmix::flexmix()] using
#'        [stats::update.formula()]. Only the left-hand side (response)
#'        of the formula is used. Default is to
#'        use the original model formula specified in [flexmix::flexmix()].
#' @param G Number of components in the mixture model (not used if zeta_p is given)
#' @param kappa_p Shrinkage parameter. Functions as if you added
#'                `kappa_p` observations according to the population mean to
#'                each component (hyperparameter for IG prior)
#' @param nu_p Degress of freedom (hyperparameter for IG prior)
#' @param zeta_p Scale (hyperparameter for IG prior). If not given the empirical
#'             variance divided by the square of the number of components
#'             is used as per Fraley and Raftery (2007).
#' @importFrom methods new
#' @importFrom mvtnorm dmvnorm
#' @import flexmix
#' @export
#' @return an object of class `"FLXC"`
#' @references
#' - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).
#'   *Ordinal Clustering with the flex-Scheme.*
#'   Austrian Journal of Statistics. _Submitted manuscript_.
#' - Fraley, C, Raftery, AE (2007)
#'   *Bayesian Regularization for Normal Mixture Estimation and Model-Based Clustering.*
#'   Journal of Classification, 24(2), 155-181
#' @example examples/regnorm.R
FLXMCregnorm <- function(formula=.~., zeta_p=NULL, kappa_p=0.01, nu_p=3, G=NULL) {
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregnorm")


    if(is.null(zeta_p) && is.null(G)) {
        stop("either parameter zeta_p or G is needed")
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

        if(is.null(zeta_p) && !is.null(G)) {
            zeta_p2 = component$var_data / G^2
        } else {
            zeta_p2 = zeta_p
        }


        nk = sum(w)
        ykbar = colSums(w*y)/nk

        muhat1 = (nk*ykbar + kappa_p*component$mu_p)/(kappa_p + nk)

        s2hat_numer1 = zeta_p2 + (kappa_p*nk)/(kappa_p+nk)*(ykbar - component$mu_p)^2
        s2hat_numer2 = rowSums(w * (t(y) - ykbar)^2)
        s2hat_denom = nu_p + nk + 3
        s2hat = (s2hat_numer1+s2hat_numer2) / s2hat_denom


        para = list(center = muhat1, s2 = s2hat, df = 2*ncol(y))
        z@defineComponent(c(para, component))
    }

    z
}


