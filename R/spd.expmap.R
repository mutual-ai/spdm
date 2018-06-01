#' Exponential map onto the space of SPD matrices
#'
#' Function projects one or more vectors onto
#' the space of SPD matrices space at a point \code{p}
#'
#' @param x Either an symmetric matrix, or a list of symmetric matrices
#' @param p An spd matrix from whose tangent space to perform the projection. If
#' unspecified, the identity matrix is chosen.
#' @details Function uses the exponential map to project a set of symmetric
#' matrices onto the space of SPD matrices at a matrix \code{p}.


spd.expmap <- function(x, p = NULL){

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
        return( p.sqrt %*% expm(p.inv.sqrt %*% x %*% p.inv.sqrt) %*% p.sqrt )
    }

    if (is.list(x)){
        ret <- lapply(x, function(i) {
            p.sqrt %*% expm(p.inv.sqrt %*% i %*% p.inv.sqrt) %*% p.sqrt
        })
        return(ret)
    }

}
