#' Internal helper functions for ctmva
#'
#' @name ctmva-internal
#' @keywords internal
NULL

#' @describeIn ctmva-internal Modified version of eigen
mygen <- function(mtrx, Im.tol=1e-15) {
    eig <- eigen(mtrx)
    evals <- eig$values
    if (max(abs(Im(evals))) > Im.tol) stop("Complex eigenvalues") else evals <- Re(evals) 
    descord <- order(evals, decreasing=TRUE)
    list(values=evals[descord], vectors=eig$vectors[ , descord])
}

#' @describeIn ctmva-internal Trapezoid weights
trapwts <- function(tt) {
    if (any(tt != sort(tt))) stop("tt should be sorted")
    nt <- length(tt)
    if (nt<3) stop("tt must have length 3 or more")
    wt <- numeric(nt)
    wt[1] <- tt[2]-tt[1]
    wt[2:(nt-1)] <- tt[-(1:2)] - tt[1:(nt-2)]
    wt[nt] <- tt[nt] - tt[nt-1]
    wt <- wt / (2 * diff(range(tt)))
    wt
}
