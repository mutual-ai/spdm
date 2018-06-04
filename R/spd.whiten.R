#' Whiten a symmetric positive definite matrix
#'
#' Approximately whitens an spd matrix x by a matrix p. Specifically,
#' conjugates x by inverse.sqrt(p), which, assuming that p is relatively close
#' to x, is close to the identity matrix.
#'
#' @param x,p Symmetric, positive-definite matrices
#' @param dewhiten If TRUE, reverse the whitening transform

spd.whiten <- function(x, p, dewhiten = F) {

    # Check input
    if (!is.spd(x) | !is.spd(p)){
        stop('x must be a symmetric matrix')
    }

    # Whiten
    if (!dewhiten){
        p.inv.sqrt <- solve(sqrt(p))
        return( p.inv.sqrt %*% x %*% p.inv.sqrt )
    } else {
        p.sqrt <- sqrt(p)
        return( p.sqrt %*% x %*% p.sqrt )
    }

}
