#' Compute means of basis functions
#'
#' Given a basis object as defined in the \pkg{fda} package (see \code{\link[fda]{basisfd}}),
#' this function simply computes the vector of means of the basis functions. Used internally.
#'
#' @param basis a basis object of class \code{"\link[fda]{basisfd}"}
#' @param rng time range. By default, the entire interval spanned by the basis. Must be left NULL for Fourier bases.
#'
#' @return  Vector of means of the basis functions
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @examples
#'
#'
#' require(fda)
#' bbasis6 <- create.bspline.basis(nbasis=6)
#' meanbasis(bbasis6)
#' meanbasis(bbasis6, c(.3,.6))
#' fbasis11 <- create.fourier.basis(nbasis=11)
#' meanbasis(fbasis11)
#'
#' @export meanbasis
meanbasis <- function(basis, rng=NULL) {
    if (basis$type %in% c("bspline", "polygonal")) {
        if (is.null(rng)) rng <- basis$rangeval
        if (rng[1]>=rng[2]) stop("'rng' must be NULL or a 2-element vector with rng[1]<rng[2]")
        inner.pts <- sort(unique(basis$params))
        inner.pts <- inner.pts[inner.pts>rng[1] & inner.pts<rng[2]]
        pts <- sort(unique(c(rng, inner.pts)))  # endpoints and internal knots
        npts <- length(pts)
        midpts <- cbind(pts[-1],pts[-npts]) %*% rbind(1:5, 5:1) / 6
        allpts0 <- c(pts, as.vector(midpts))
        if (npts>2) allpts0 <- c(allpts0, pts[2:(npts-1)])  # add a second copy of each internal knot
        allpts <- sort(allpts0)
        mult <- rep(diff(pts), each=7) * rep(c(41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140),(npts-1)) / 6
        dmult <- diag(mult)
        return(colSums(dmult %*% fda::eval.basis(allpts, basis)) / diff(rng))
    }
    else if (basis$type=="fourier") {
        if (!is.null(rng)) stop("Argument 'rng' implemented only for B-spline bases")
        return(c(1/sqrt(diff(basis$rangeval)),rep(0,(basis$nbasis-1))))
    }
    else stop("'basis' must be a B-spline or Fourier basis.")
}
