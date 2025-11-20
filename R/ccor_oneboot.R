#' Curve correlation for a single bootstrap sample
#'
#' This internal function is called repeatedly by \code{\link{ccor_boot}}. It is not ordinarily called by users.
#'
#' @param y,time,curve,k,min.overlap,min.n see \code{\link{ccor}}
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>, Noemi Foa, Dror Arbiv and Biplab Paul <paul.biplab497@gmail.com>
#' @keywords internal
#' @export

ccor_oneboot <-
function(y, time, curve=NULL, k = 15, min.overlap = 0, min.n=8) {

  boot_indices <- sort(sample.int(length(time), replace = TRUE))
  tboot <- time[boot_indices]
  yboot <- if (is.vector(y)) y[boot_indices] else y[boot_indices, ]
  if (!is.null(curve)) curveboot <- curve[boot_indices]

  ntime <- length(unique(tboot))
  if (k > ntime) {
    k <- min(k, ntime)
    warning(paste("k truncated to", ntime))
  }

  ccor(y=yboot, time=tboot, curve=curveboot, k = k, min.overlap = min.overlap, min.n=min.n)$estcor
}
