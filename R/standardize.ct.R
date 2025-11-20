#' Center and scale a continuous-time multivariate data set
#'
#' Subtracts the (continuous-time) mean and divides by the (continuous-time) standard
#' deviation of each of the variables. This is the continuous-time analogue
#' of taking an \eqn{n \times p} data matrix, subtracting the mean of each column, and
#' dividing by the standard deviation of each column, as is done by
#' \code{\link{scale}(\dots, center=TRUE, scale=TRUE)}.
#'
#' @param fdobj continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @return A standardized (centered and scaled) version of the input data.
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il> and Biplab Paul <paul.biplab497@gmail.com>
#' @seealso  \code{\link{center.ct}}
#' @export standardize.ct
standardize.ct<- function(fdobj) {
    cfd <- center.ct(fdobj)
    covmat <- cov.ct(fdobj)
    fd(coef = cfd$coefs %*% diag(1/sqrt(diag(covmat))),
       basisobj = fdobj$basis, fdnames=fdobj$fdnames)
}
