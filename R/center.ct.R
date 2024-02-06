#' Center a continuous-time multivariate data set
#' 
#' Subtracts the (continuous-time) mean of each of the variables. This is analogous
#' to column-centering an \eqn{n \times p} data matrix. 
#' 
#' @param fdobj continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @return A centered version of the input data.
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#' @seealso  \code{\link{standardize.ct}} 
#' @export center.ct
center.ct <- function(fdobj) {
    phibar <- meanbasis(fdobj$basis)
    nbasis <- fdobj$basis$nbasis
    if (fdobj$basis$type %in% c("bspline", "polygonal")) {
        w <- rep(1,nbasis)
    }
    else if (fdobj$basis$type=="fourier") {
        period <- diff(fdobj$basis$rangeval)
        w <- c(sqrt(period),rep(0,(nbasis-1)))
    }
    else stop("The basis for 'fdobj' must be a B-spline or Fourier basis.")
    newcoefmat <- fdobj$coefs - w %*% crossprod(phibar, fdobj$coefs)
    fd(coef = newcoefmat, basisobj = fdobj$basis, fdnames=fdobj$fdnames)
}