#' FlexMix Driver for Regularized Multinomial Mixtures
#'
#' This model driver can be used to cluster data using a multinomial
#' distribution.
#' 
#' Using a regularization parameter `alpha2` greater than zero
#' acts as adding `alpha2` observations conforming to the population
#' mean to each component. This can be used to avoid degenerate
#' solutions. It also has the effect
#' that clusters become more similar to each other the larger
#' `alpha2` is chosen. For small values it is mostly negligible however.
#'
#' For regularization we compute the MAP estimates for the multinomial
#' distribution using the Dirichlet distribution as prior, which is 
#' the conjugate prior. The parameters of this prior are selected to 
#' correspond to the marginal distribution of the variable across all
#' observations.
#'
#' @param size values are assumed to be integers in `1:size`
#' @param alpha2 Regularization parameter. Can be regarded the same as
#'  adding `alpha2` observations conforming to the population mean to each
#'  component.
#' @param formula A formula which is interpreted relative to the formula
#'        specified in the call to `flexmix` using `update.formula`. Only the
#'        left-hand side (response) of the formula is used. Default is to
#'        use the original `flexmix` model formula.
#' @return an object of class FLXC
#' @export
#' @references
#' - Galindo Garre, F, Vermunt, JK (2006).
#'   *Avoiding Boundary Estimates in Latent Class Analysis by Bayesian Posterior Mode Estimation*
#'   Behaviormetrika, 33, 43-59.
#' - Ernst, D, Ortega Menjivar, L, Scharl T, Grün, B (2025).
#'   *Ordinal clustering with the flex-Scheme.*
#'   Austrian Statistics Journal. _Submitted manuscript_.
#' @export
#' @example examples/multinom.R
FLXMCregmultinom = function(formula=.~., size, alpha2=0) {
    stopifnot(length(size) == 1)
    stopifnot(length(alpha2) == 1)

    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregmultinom")

    .as01 = function(y) {
        # as matrix in 0/1 coding
        yd = lapply(seq_len(ncol(y)), \(col) {
            #uvalues = sort(unique(y[, col]))
            uvalues = seq_len(size)
            cns = paste0(colnames(y)[col], uvalues)

            xx = lapply(uvalues, \(value) {
                as.integer(y[, col] == value)
            }) |> do.call(cbind, args=_) |>
               `colnames<-`(cns)
            xx
        }) |> do.call(cbind, args=_)
        yd
    }

    z@preproc.y <- function(y) {
        if (any(y < 1, na.rm=TRUE)) {
            stop("values < 1 not allowed (values need to be in 1:size)")
        }

        if(any(y > size)) {
            stop("values larger than size not allowed (values need to be in 1:size)")
        }

        y
    }

    defineComponent <- function(component) {
        probs = component$probs
        predict <- function(x, ...) {
            stop("not implemented")
        }

        logLik <- function(x, y) {
            probmat = lapply(seq.int(ncol(y)), \(col) {
                probs[, col][y[, col]]
            }) |> do.call(cbind, args=_) |>
                `colnames<-`(colnames(y))
            rowSums(log(probmat))
        }

        new("FLXcomponent",
            parameters=component,
            logLik=logLik, df=component$df, 
            predict=predict)
    }

    z@fit <- function(x, y, w, component) {
        if(length(component) == 0) {
            component$yd = .as01(y)
            component$ymarg = colMeans(component$yd) |> unname()
            component$b_alpha = component$ymarg*alpha2
            component$b_beta = (1-component$ymarg)*alpha2
        }

        p = with(component, (b_alpha + colSums(w*yd, na.rm=TRUE)) / (b_alpha+b_beta+sum(w))) |>
            matrix(ncol=ncol(y), nrow=size)
        component$probs = p
        component$df = ncol(y)
        defineComponent(component)
    }

    z
}
