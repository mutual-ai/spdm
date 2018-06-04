#' Exponential map onto the space of SPD matrices
#'
#' Function projects one or more vectors onto
#' the space of SPD matrices space at a point \code{p}
#'
#' @param x A symmetric matrix
#' @param p An spd matrix from whose tangent space to perform the projection. If
#' unspecified, the identity matrix is chosen.
#' @details Function uses the exponential map to project a symmetric
#' matrix onto the space of SPD matrices at a matrix \code{p}.


spd.expmap <- function(x, p = NULL){

    # Check x input
    if (!'s.mat' %in% input.type(x)){
        stop('x must be a symmetric matrix')
    }

    # Set and check p
    if (is.null(p)){
        p <- diag(rep(1, dim(x)[1]))
    } else if (!'spd.mat' %in% input.type(p)){
        stop('p must be a symmetric positive definite matrix')
    }

    p.sqrt     <- sqrtm(p)
    p.inv.sqrt <- solve(p.sqrt)
    return(p.sqrt %*% expm(p.inv.sqrt %*% x %*% p.inv.sqrt) %*% p.sqrt)

}
