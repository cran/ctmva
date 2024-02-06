#' Continuous-time covariance or cross-covariance matrix
#'
#' Computes the covariance matrix of a continuous-time multivariate
#' data set represented as an \code{\link[fda]{fd}} object; or the
#' cross-covariance matrix of two such data sets.
#'
#' @param fdobj1  continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param fdobj2  an optional second data set
#' @param common_trend  logical: centering with respect to the mean function if \code{TRUE},
#' without centering if \code{FALSE} (the default)
#' @return  A matrix of (cross-) covariances
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso \code{\link{cor.ct}}
#' @examples
#'
#'# see example for cor.ct, which works similarly
#'
#' @export cov.ct
cov.ct <- function(fdobj1, fdobj2=fdobj1, common_trend=FALSE) {
    if (common_trend) {
      fdobj1 <- fda::center.fd(fdobj1)
      fdobj2 <- fda::center.fd(fdobj2)
    }
    rng1 <- fdobj1$basis$rangeval; rng2 <- fdobj2$basis$rangeval
    overlap <- c(max(rng1[1],rng2[1]), min(rng1[2],rng2[2]))
    P0 <- inprod.cent(fdobj1$basis, fdobj2$basis)  # range of integration automatically set to 'overlap'
    Phi1 <- t(fdobj1$coef)
    Phi2 <- t(fdobj2$coef)
    matt <-  Phi1 %*% P0 %*% t(Phi2) / diff(overlap)

    if(!is.null(fdobj1$fdnames$reps)) rownames(matt) <- fdobj1$fdnames$reps
    else rownames(matt) <- paste0("Comp.", 1:ncol(fdobj1$coef))

    if(!is.null(fdobj2$fdnames$reps)) colnames(matt) <- fdobj2$fdnames$reps
    else colnames(matt) <- paste0("Comp.", 1:ncol(fdobj2$coef))

    if (length(matt)==1) matt <- matt[1,1]

    return(matt)
}
