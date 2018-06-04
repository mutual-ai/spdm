#' Log-map onto the tangent space
#'
#' Function projects one or more symmetric, positive-definite matrices onto
#' the tangent space at a point \code{p}
#'
#' @param x An spd matrix
#' @param p An spd matrix on whose tangent space to perform the projection. If
#' unspecified, the identity matrix is chosen.
#' @details Function uses the log-map to project a symmetric, positive-definite
#' matrix onto the tangent space of the positive cone at a matrix \code{p}.


spd.logmap <- function(x, p = NULL){

    # Check x input
    if (!'spd.mat' %in% input.type(x)){
        stop('x must be a symmetric positive definite matrix')
    }

    # Set and check p
    if (is.null(p)){
        p <- diag(rep(1, dim(x)[1]))
    } else if (!'s.mat' %in% input.type(p)){
        stop('p must be a symmetric positive definite matrix')
    }

    p.sqrt     <- sqrtm(p)
    p.inv.sqrt <- solve(p.sqrt)
    return( p.sqrt %*% logm(p.inv.sqrt %*% x %*% p.inv.sqrt) %*% p.sqrt )

}
