#' Log-map onto the tangent space
#'
#' Function projects one or more symmetric, positive-definite matrices onto
#' the tangent space at a point \code{p}
#'
#' @param x Either an spd matrix, or a list of spd matrices
#' @param p An spd matrix on whose tangent space to perform the projection. If
#' unspecified, the identity matrix is chosen.
#' @details Function uses the log-map to project a set of symmetric, positive-definite
#' matrices onto the tangent space of the positive cone at a matrix \code{p}.


spd.logmap <- function(x, p = NULL){

    # Verify input
    input.check <- input.type(x)
    if (!is.null(p)){
        p.check <- input.type(p)
        if (!'spd.mat' %in% p.check) {
            stop('p is not spd')
        }
    }

    if (!'spd.mat' %in% input.type(x) | !'spd.mat' %in% input.type(y)){
        return('Both inputs must be positive definite matrices')
    }

    if (is.null(p)){
        if (is.matrix(x)){
            p <- diag(rep(1, dim(x)[1]))
        } else if (is.list(x)) {
            p <- diag(rep(1, dim(x[[1]])[1]))
        }
    }

    p.sqrt     <- sqrtm(p)
    p.inv.sqrt <- solve(p.sqrt)

    if (is.matrix(x)){
        return( p.sqrt %*% logm(p.inv.sqrt %*% x %*% p.inv.sqrt) %*% p.sqrt )
    }

    if (is.list(x)){
        ret <- lapply(x, function(i) {
            p.sqrt %*% logm(p.inv.sqrt %*% i %*% p.inv.sqrt) %*% p.sqrt
        })
        return(ret)
    }

}
