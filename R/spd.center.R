#' Center a set of spd matrices
#'
#' Functions centers a set of symmetric, positive-definite matrices so that they
#' have mean equal to the identity matrix.
#'
#' @param x Either a list of symmetric, positive-definite matrices, an individual amtrix
#' @param y Value(s) for centering. Either a list of spd matrices of the same
#' length as x, a single spd matrix, or NULL. If NULL, uses the mean of \code{x}
#' @details Functions transports the matrices in \code{x} so that \code{y} lies at
#' the identity. If \code{x} is a list and \code{y = NULL} (default), function sets
#' \code{y} equal to the mean of \code{x}. Otherwise, \code{y} can be set to any
#' spd matrix. When \code{x} is a matrix, \code{y} must be provided. If \code{y}
#' is a list of the same length as {x}, the function centers each element of
#' \code{x} by the corresponding element of \code{y}.
#'
#' Centering is accomplished by the parallel transport specified by XXX (2XXX), by
#' first computing the inverse square root of \code{y}, and then centering \code{x}
#' by \code{y.inv.sqrt x y.inv.sqrt}.
#'
#' @return An object of the same type as \code{x}, containing the centered matrices..

spd.center <- function(x, y = NULL){

    if (is.list(x)){

        # If y is a list, do elementwise
        if (is.list(y)){
            if (length(x) != length(y)){
                stop('If a list, y must have the same length as x')
            }

            # inv-sqrm
            y.inv.sqrt <- lapply(y, function(i) sqrtm(solve(i)))

            # Center
            x.centered <- lapply(1:length(x), function(i){
                y.inv.sqrt[i] %*% x[i] %*% y.inv.sqrt[i]
            })

        }

        # If no center specified, use mean
        if (is.null(y)){
            y <- spd.mean(x)
        }

        # inv-sqrm
        y.inv.sqrt <- sqrtm(solve(y))

        # Center
        x.centered <- lapply(x, function(i)
            y.inv.sqrt %*% i %*% y.inv.sqrt)

    } else if (is.matrix(x)) {

        # If no center specified, abort
        if (is.null(y)){
            stop('Must provide y if x is a single matrix')
        }

        # inv-sqrm
        y.inv.sqrt <- sqrtm(solve(y))

        # Center
        x.centered <- y.inv.sqrt %*% x %*% y.inv.sqrt
    }

    return(x.centered)
}
