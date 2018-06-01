#' Compute the mean of a set of spd matrices
#'
#' Function computes the mean of a set of symmetric positive-definite matrices.
#' Several methods are implemented, as described in \code{Details}.
#'
#' @param x A list of symmetric, positive-definite matrices
#' @param method The type of mean to compute. See details
#' @param ... Further arguments. See details
#' @details Function computes the mean of a set of symmetrix, positive-definite
#' matrices. Several methods are implemented:
#' \itemize{
#'  \item{"euclidean": }{The ordinary arithmetic mean -- the sum of the matrices in x,
#'  divided by the number of matrices. This is guaranteed to be an spd matrix, but
#'  does not necessarily preserve the spectral characteristics of the individual matrices.}
#'  \item{"logeuclidean": }{Computed by taking the arithmetic mean of the logarithms
#'  of the matrices in \code{x}, and then projecting back onto the space of spd matrices.
#'  In general, better behaved than the arithmetic mean.}
#'  \item{"riemannian": }{The geometric mean proposed by Barachant et al. (2013).
#'  Is approximated iteratively, and accepts two additional argument: \code{max.iter},
#'  giving the maximum number of iterations (defaults to 20), and \code{tol}, giving
#'  the maximum Frobenius norm of the approximation error (defaults to 0.1)}
#' }
#' @return The mean of the matrices in \code{x}.

spd.mean <- function(x, method = 'euclidean'){

    if (method == 'euclidean'){
        return(Reduce(`+`, x)/length(x))
    }

    if (method == 'logeuclidean'){
        logmean <- Reduce(`+`, spd.logmap(x))/length(x)
        return(expm(logmean))
    }

    if (method == 'riemannian'){

        # Initialize to arithmetic mean
        init <- Reduce(`+`, x) / length(x)

        # Iterate until convergence
        err  <- 2*tol
        iter <- 1
        while (err > tol & iter <= max.iter){

            # Project onto tangent space at current mean estimate
            proj <- spd.logmap(x, p = init)
            tan.mean <- Reduce(`+`, proj) / length(proj)

            # Update mean estimate
            init <- spd.expmap(tan.mean, p = init)

            # Calculate error
            err <- norm(tan.mean, type = 'F')
            iter <- iter + 1
        }

        return(init)

    }
}
