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
#' @param r Number of different categories. Values are assumed to be
#'     integers in `1:r`. Default `NULL` implies that the number of
#'     different categories is inferred columnwise by the maximum
#'     value observed.
#' @param alpha A non-negative scalar acting as regularization
#'     parameter. Can be regarded as adding `alpha` observations equal
#'     to the population mean to each component.
#' @param formula A formula which is interpreted relative to the
#'     formula specified in the call to [flexmix::flexmix()] using
#'     [stats::update.formula()]. Only the left-hand side (response)
#'     of the formula is used. Default is to use the original model
#'     formula specified in [flexmix::flexmix()].
#' @return an object of class `"FLXC"`
#' @export
#' @references - Galindo Garre, F, Vermunt, JK (2006).  *Avoiding
#'     Boundary Estimates in Latent Class Analysis by Bayesian
#'     Posterior Mode Estimation* Behaviormetrika, 33, 43-59.  -
#'     Ernst, D, Ortega Menjivar, L, Scharl, T, Grün, B (2025).
#'     *Ordinal Clustering with the flex-Scheme.* Austrian Journal of
#'     Statistics. _Submitted manuscript_.
#' @export
#' @importFrom utils head
#' @example examples/multinom.R
FLXMCregmultinom = function(formula=.~., r = NULL, alpha=0) {
    stopifnot(is.numeric(alpha), length(alpha) == 1, alpha >= 0)
    b_alpha <- 0

    z <- new("FLXMC", weighted=TRUE, formula=formula,
             name="FLXMCregmultinom")

    .as01 <- function(y, r) {
        # as matrix in 0/1 coding
        yd = lapply(seq_len(ncol(y)), \(col) {
            uvalues = seq_len(r[col])
            cns = paste0(colnames(y)[col], uvalues)

            xx = lapply(uvalues, \(value) {
                as.integer(y[, col] == value)
            }) |> do.call(cbind, args=_) |>
               `colnames<-`(cns)
            xx
        }) |> do.call(cbind, args=_)
        attr(yd, "y") <- y
        attr(yd, "ymarg") <- colMeans(yd) |> unname()

        csr <- cumsum(r)
        csrlag <- c(1, head(csr, -1)+1)        
        attr(yd, "csr") <- csr
        attr(yd, "csrlag") <- csrlag
        
        yd
    }

    z@preproc.y <- function(y) {
        if (anyNA(y)) 
            stop("NAs are not allowed")
        if (any(y < 1)) {
            stop("values < 1 not allowed (values need to be in 1:r)")
        }
        if (is.null(r)) {
            r <- apply(y, 2, max)
        } else {
            r <- as.integer(r)
            r <- rep(r, length = ncol(y))
            if (any(apply(y, 2, max) > r)) {
                stop("values larger than r not allowed (values need to be in 1:r)")
            }
        }
        stopifnot(r >= 2)

        .as01(y, r)
    }

    z@defineComponent <- function(para) {
        predict <- function(x, ...) {
            stop("not implemented")
        }
        
        logLik <- function(x, y) {
            y <- attr(y, "y")
            probmat <- lapply(seq_along(para$prob), \(col) {
                log(para$prob[[col]])[y[, col]]
            }) |> do.call(cbind, args=_) |>
                `colnames<-`(colnames(y))
            rowSums(probmat)
        }

        new("FLXcomponent",
            parameters=list(p = para$prob),
            logLik=logLik, df=para$df, 
            predict=predict)
    }

    z@fit <- function(x, y, w) {
        if (alpha > 0) {
            b_alpha <- attr(y, "ymarg") * alpha
        } 

        p <- (b_alpha + colSums(w*y)) / (alpha + sum(w))
        pl <- mapply(\(a,b) {
            p[a:b]
        }, attr(y, "csrlag"), attr(y, "csr"), SIMPLIFY = FALSE)

        para <- list(prob = pl, df = sum(lengths(pl) - 1))
        z@defineComponent(para)
    }

    z
}
