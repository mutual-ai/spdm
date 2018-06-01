#' Distance between two symmetric, positive-definite matrices
#'
#' Function implements several distance measures between symmetric,
#' positive-definite matrices, not all of which are proper metrics.
#'
#' @param x,y Symmetric, positive-definite matrices
#' @param method The distance measure. See details.
#' @details Allowable distance measures are
#' \itemize{
#'  \item{"frobenius": }{The Frobenius norm of the difference \code{x-y}. Not affinely invariant.}
#'  \item{"cholesky": }{The Frobenius norm of the difference between the cholesky factors
#'  of \code{x} and \code{y}. Not affinely invariant.}
#'  \item{"affine": }{An affinely invariant distance measure, computed by taking
#'  the Frobenius norm of the projection of \code{y} into the tangent space at \code{x}}
#'  \item{"riemannian": }{The Riemmanian distance proposed by Barachant, et al. (2013)}
#' }

spd.dist <- function(x, y, method = 'riemannian'){

    if (!'spd.mat' %in% input.type(x) | !'spd.mat' %in% input.type(y)){
        return('Both inputs must be positive definite matrices')
    }

    if (method == 'frobenius'){
        return(norm(x - y, type = 'F'))
    }

    if (method == 'cholesky'){
        return(norm(chol(x) - chol(y), type = 'F'))
    }

    if (method == 'affine'){
        isr <- solve(sqrtm(x))
        return(norm(logm(isr %*% y %*% isr), type = 'F'))
    }

    if (method == 'riemannian'){
        return(norm(logm(solve(x) %*% y), type = 'F'))
    }
}
