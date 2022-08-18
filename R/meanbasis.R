#' Compute means of basis functions
#'
#' Given a basis object as defined in the \pkg{fda} package (see \code{\link[fda]{basisfd}}),
#' this function simply computes the vector of means of the basis functions. Used internally.
#'
#' @param basis a basis object of class \code{"\link[fda]{basisfd}"}
#' @return  Vector of means of the basis functions
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @examples
#'
#'
#' require(fda)
#' fbasis11 <- create.fourier.basis(nbasis=11)
#' zapsmall(meanbasis(fbasis11))    # the sine functions have mean 0
#'
#' @export meanbasis
meanbasis <- function(basis) {
    pts <- sort(unique(c(basis$rangeval, basis$params))); npts <- length(pts)
    midpts <- cbind(pts[-1],pts[-npts]) %*% rbind(1:5, 5:1) / 6
    allpts0 <- c(pts, as.vector(midpts))
    if (npts>2) allpts0 <- c(allpts0, pts[2:(npts-1)])
    allpts <- sort(allpts0)
    mult <- rep(diff(pts), each=7) * rep(c(41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140),(npts-1)) / 6
    dmult <- diag(mult)
    colSums(dmult %*% eval.basis(allpts, basis)) / diff(basis$rangeval)
}
