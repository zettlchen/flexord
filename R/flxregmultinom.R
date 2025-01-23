#' @export
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
