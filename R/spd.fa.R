#' Fractional anisotropy of a symmetrix, positive-definite matrix
#'
#' @param x A symmetric, positive-definite matrix

spd.fa <- function(x){

    # Check x input
    if (!'s.mat' %in% input.type(x)){
        stop('x must be a symmetric matrix')
    }

    k <- dim(x)[1]
    l <- eigen(x, symmetric = T, only.values = T)$values
    l.bar <- mean(l)
    fa <- sqrt( (k)/(k-1) *   sum( (l - l.bar)^2) / sum(l^2)  )
    return(fa)

}
