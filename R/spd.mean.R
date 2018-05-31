#' Compute the mean of a set of spd matrices
#'
#' Function computes the mean of a set of symmetric positive-definite matrices
#' by taking the logarithm of each amtrix, taking the mean, and then transforming
#' back using the matrix exponential. By default, the function tries to diagonalize
#' each argument, and then uses the method by Higham (2008) if that fails.
#'
#' @param x A list of symmetric, positive-definite matrices
#' @details The function uses the \code{logm} and \code{expm} functions in the
#' \code{expm} package. It first attempts to diagonalize each matrix using
#' \code{method = 'Eigen'}, switching to \code{meth0d = 'Higham08'} if that fails.
#'
#' @return The mean of the matrices in \code{x}.

spd.mean <- function(x){

    # Matrix logarithms
    lx <- lapply(x, function(i) {
        tryCatch(logm(i, method = 'Eigen'),
                 error = function(e) {
                     logm(i)
                 })
    })

    # Mean and exponential
    lm <- Reduce(`+`, lx) / length(lx)
    m  <- tryCatch(expm(lm, method = 'Eigen'),
                   error = function(e) expm(lm))

    return(m)
}
