#' Curve correlation
#'
#' Inputs raw data representing two curves, applies penalized B-spline
#' smoothing to the two curves, and computes the curve correlation between
#' them via a call to \code{\link{cor.ct}}.
#'
#' If \code{y} is a two-column matrix, the two curves are observed at the time points given by
#' \code{time}; in this case \code{length(time)} must equal \code{nrow(y)}, and \code{curve} is
#' ignored. If \code{y} is a vector, it must have the same length as both \code{time} and \code{curve}.
#' In this case \code{y} contains the observations on both curves, while elements of \code{time} and \code{curve}
#' identify the observation time and the curve being observed, respectively.
#'
#' @param y either a vector or a two-column matrix, without missing values; see Details
#' @param time a vector of time points
#' @param curve curve indicator; see Details
#' @param k number of B-spline basis functions
#' @param min.overlap minimum overlap of the two curves' time ranges
#' @param min.n minimum number of observations per curve
#' @return A list with components
#' \item{y,time }{the supplied \code{y} and \code{time}}
#' \item{mod1,mod2 }{models for the two curves, outputted by \code{\link[mgcv]{gam}}}
#' \item{fd1,fd2 }{functional data objects (see \code{\link[fda]{fd}}) for the two curves}
#' \item{estcor }{estimated curve correlation}
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>, Noemi Foa, Dror Arbiv and Biplab Paul <paul.biplab497@gmail.com>
#' @seealso \code{\link{cor.ct}}, \code{\link[mgcv]{b.spline}}
#' @examples
#'
#' # see example for ccor_posim
#'
#' @importFrom MASS ginv
#' @importFrom stats fitted
#' @importFrom fda create.bspline.basis
#' @export

ccor <-
function(y, time, curve=NULL, k = 15, min.overlap = 0, min.n=8) {
  # Check that time is a vector
  if (!is.vector(time)) stop("time must be a vector")

  # Check that y is either a vector or a 2-column matrix
  if (!is.vector(y) && !(is.matrix(y) && ncol(y) == 2)) {
    stop("y must be either a vector or a 2-column matrix")
  }

  # If y is a vector, check the length of time and curve
  if (is.vector(y)) {
    if (length(time) != length(y) || length(curve) != length(y)) {
      stop("If y is a vector, time and curve must be vectors of the same length as y")
    }
    curve <- factor(curve)
    if (length(levels(curve)) != 2)
        stop("If y is a vector, curve must take exactly two values")
    ylist <- list(y[curve == levels(curve)[1]], y[curve == levels(curve)[2]])
  }

  # If y is a 2-column matrix
  if (is.matrix(y)) {
    if (nrow(y) != length(time)) stop("If y is a matrix, its number of rows must match the length of time")
    ylist <- list(y[, 1], y[, 2])
  }

  n1 <- sum(!is.na(ylist[[1]]))
  n2 <- sum(!is.na(ylist[[2]]))
  if (min(n1, n2) < min.n) stop("Not enough data")

  tt <- if (is.vector(y)) time[curve == levels(curve)[1]] else time
  mod1 <- mgcv::gam(ylist[[1]] ~ s(tt, bs="bs", k=k), method = "REML")
  wots1 <-  mod1$smooth[[1]]$knots
  nots1 <- sort(c(range(tt), wots1[wots1>min(tt) & wots1<max(tt)]))
  bsb1 <- fda::create.bspline.basis(range(tt), breaks=nots1)
  ginvB1 <- ginv(fda::eval.basis(tt, bsb1))
  nucoef1 <- ginvB1 %*% fitted(mod1)
  est_y1 <- fd(coef = nucoef1, basisobj = bsb1)

  if (is.vector(y)) tt <- time[curve == levels(curve)[2]]
  mod2 <- mgcv::gam(ylist[[2]] ~ s(tt, bs="bs", k=k), method = "REML")
  wots2 <-  mod2$smooth[[1]]$knots
  nots2 <- sort(c(range(tt), wots2[wots2>min(tt) & wots2<max(tt)]))
  bsb2 <- create.bspline.basis(range(tt), breaks=nots2)
  ginvB2 <- ginv(fda::eval.basis(tt, bsb2))
  nucoef2 <- ginvB2 %*% fitted(mod2)
  est_y2 <- fd(coef = nucoef2, basisobj = bsb2)

  # Check overlap in time ranges
  if (min(bsb1$rangeval[2], bsb2$rangeval[2]) - max(bsb1$rangeval[1], bsb2$rangeval[1]) < min.overlap)
    stop("Insufficient overlap of time ranges")

  res <- list(y = y, time = time,
              mod1 = mod1, mod2 = mod2,
              fd1 = est_y1, fd2 = est_y2,
              estcor = cor.ct(est_y1, est_y2))

  res
}
