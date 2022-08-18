#'  Centered inner product matrix for a basis or pair of bases
#'
#' Most methods of continous-time multivariate analysis require a matrix of
#' inner products of pairs of functions from a basis, such as a B-spline basis,
#' or pairs consisting of one function from each of two bases. This function computes
#' such matrices via 7-point Newton-Cotes integration, which is exact for cubic
#' B-splines.
#'
#' @param basis1  basis object from the \code{\link[fda]{fda}} package.
#' @param basis2  an optional second basis
#' @param rng  range (of times) spanned by the basis
#' @return  Matrix of inner products of each pair of centered basis functions.
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso  \code{\link[fda]{create.bspline.basis}} from package \code{\link[fda]{fda}}, for the most commonly used basis object type.
#' @examples
#'
#'
#' require(fda)
#' basis8 <- create.bspline.basis(nbasis=8)
#' inprod.cent(basis8)
#'
#' @export inprod.cent
inprod.cent <-
function(basis1, basis2=basis1, rng=NULL) {
    if (is.null(rng)) rng <- c(max(basis1$rangeval[1],basis2$rangeval[1]), min(basis1$rangeval[2],basis2$rangeval[2]))
    if (rng[1]>=rng[2]) stop("'rng' must be a 2-element vector with rng[1]<rng[2]")
    inner.pts <- sort(unique(c(basis1$params, basis2$params)))
    inner.pts <- inner.pts[inner.pts>rng[1] & inner.pts<rng[2]]
    pts <- sort(unique(c(rng, inner.pts))); npts <- length(pts)
    midpts <- cbind(pts[-1],pts[-npts]) %*% rbind(1:5, 5:1) / 6
    allpts0 <- c(pts, as.vector(midpts))
    if (npts>2) allpts0 <- c(allpts0, pts[2:(npts-1)])
    allpts <- sort(allpts0)
    mult <- rep(diff(pts), each=7) * rep(c(41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140),(npts-1)) / 6
    dmult <- diag(mult)
    evalmat1 <- eval.basis(allpts, basis1)
    int1 <- colSums(dmult %*% evalmat1)
    evalmat2 <- eval.basis(allpts, basis2)
    int2 <- colSums(dmult %*% evalmat2)
    inprod <- t(evalmat1) %*% dmult %*% evalmat2
    res <- inprod - outer(int1, int2)/diff(rng)
    if (length(all.equal(basis1, basis2))==1) res <- (res + t(res)) / 2
    res
}
