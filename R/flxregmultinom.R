#' FlexMix Driver for Regularized Multinomial Mixtures
#'
#' This model driver can be used to cluster data using a multinomial
#' distribution.
#' 
#' Using a regularization parameter `alpha` greater than zero
#' acts as adding `alpha` observations conforming to the population
#' mean to each component. This can be used to avoid degenerate
#' solutions. It also has the effect
#' that clusters become more similar to each other the larger
#' `alpha` is chosen. For small values it is mostly negligible however.
#'
#' For regularization we compute the MAP estimates for the multinomial
#' distribution using the Dirichlet distribution as prior, which is 
#' the conjugate prior. The parameters of this prior are selected to 
#' correspond to the marginal distribution of the variable across all
#' observations.
#'
#' @param r Number of different categories. Values are assumed to be integers in `1:r`.
#' @param alpha A non-negative scalar acting as regularization
#'     parameter. Can be regarded as adding `alpha` observations
#'     equal to the population mean to each component.
#' @param formula A formula which is interpreted relative to the formula
#'        specified in the call to [flexmix::flexmix()] using
#'        [stats::update.formula()]. Only the
#'        left-hand side (response) of the formula is used. Default is to
#'        use the original model formula specified in [flexmix::flexmix()].
#' @return an object of class `"FLXC"`
#' @export
#' @references
#' - Galindo Garre, F, Vermunt, JK (2006).
#'   *Avoiding Boundary Estimates in Latent Class Analysis by Bayesian Posterior Mode Estimation*
#'   Behaviormetrika, 33, 43-59.
#' - Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).
#'   *Ordinal Clustering with the flex-Scheme.*
#'   Austrian Journal of Statistics. _Submitted manuscript_.
#' @export
#' @importFrom utils head
#' @example examples/multinom.R
FLXMCregmultinom = function(formula=.~., r, alpha=0) {
    stopifnot(is.numeric(alpha), length(alpha) == 1, alpha >= 0)

    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregmultinom")

    .as01 = function(y) {
        if(length(r) == 1) {
            r = rep(r, ncol(y))
        }

        # as matrix in 0/1 coding
        yd = lapply(seq_len(ncol(y)), \(col) {
            #uvalues = sort(unique(y[, col]))
            uvalues = seq_len(r[col])
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
            stop("values < 1 not allowed (values need to be in 1:r)")
        }

        if(any(apply(y, 2, max) > r)) {
            stop("values larger than r not allowed (values need to be in 1:r)")
        }

        if(length(r) != 1 && length(r) != ncol(y)) {
            stop("r must be either a scalar or a vector of length ncol(y)")
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
                #probs[, col][y[, col]]
                probs[[col]][y[, col]]
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
            component$b_alpha = component$ymarg*alpha
            component$b_beta = (1-component$ymarg)*alpha
        }

        p = with(component, (b_alpha + colSums(w*yd, na.rm=TRUE)) / (b_alpha+b_beta+sum(w)))

        csr = cumsum(r)
        csrlag = c(1, head(csr, -1)+1)

        # probabilities as list, so it's easier to compute the loglik
        pl = mapply(\(a,b) {
            p[a:b]
        }, csrlag, csr, SIMPLIFY=FALSE)

        component$probs = pl
        #component$df = ncol(y)
        component$df = sum(r-1)
        defineComponent(component)
    }

    z
}
