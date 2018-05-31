#' Regularized covariance estimation
#'
#' Function performs regularized covariance estimation using either soft-thresholding
#' or the graphical lasso, involving l1 penalized estimation of the inverse
#' covariance (precision) matrix. In either case, the optimal tuning parameter
#' is selected by cross-validation. Currently, the function's internals are set
#' to my own default settings, for simplicity.
#'
#' @param y A data matrix, where rows are observations and columns are variables
#' @param method Either soft thresholding ('threshold') or graphical lasso ('sparse')
#' @return A covariance matrix


spd.estimate <- function(y, method = 'threshold'){

    if (method == 'threshold'){
        est <- regular.CV(y, k.grid = seq(0,1,.1),
                          method = "SoftThresholding", fold = 5)
        fit <- soft.thresholding(cov(y), c = est$CV.k[1])
    } else if (method == 'sparse') {
        n <- dim(y)[1]
        est <- huge(y, method = 'glasso', scr = T, scr.num = n/log(n),
                    cov.output = T)
        sel <- huge.select(est, criterion = 'stars')
        fit <- as.matrix(sel$opt.cov)
    }
    return(fit)
}
