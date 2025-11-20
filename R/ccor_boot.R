#' Bootstrap confidence interval estimation for curve correlation
#'
#' Inputs raw data representing two curves, and computes a bootstrap confidence interval
#' for the curve correlation between them.
#'
#' @param y,time,curve,k,min.overlap,min.n see \code{\link{ccor}}
#' @param ndraw number of bootstrap samples
#' @param conf confidence level
#' @return A list with components
#' \item{cor }{curve correlation (for the full data)}
#'  \item{boot.cor }{curve correlations for the \code{ndraw} bootstrap samples}
#'  \item{ci }{confidence interval for the curve correlation}
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>, Noemi Foa, Dror Arbiv and Biplab Paul <paul.biplab497@gmail.com>
#' @seealso \code{\link{ccor}}, \code{\link{ccor_posim}}
#' @examples
#'
#' # see example for ccor_posim
#'
#' @export

ccor_boot <-
function(y, time, curve=NULL, k = 15, ndraw = 299, conf = 0.95, min.overlap = 0, min.n=8) {

  all_estcor <- numeric(ndraw)

  cor_fulldat <- ccor(y = y, time = time, curve=curve, k = k, min.overlap = min.overlap, min.n=min.n)

  for (nbt in 1:ndraw)
    all_estcor[nbt] <- ccor_oneboot(y = y, time = time, curve=curve, k = k, min.overlap = min.overlap, min.n=min.n)

  return(list(cor = cor_fulldat$estcor,
              boot.cor = all_estcor,
              ci = sort(all_estcor)[c(floor((ndraw + 1) * (1 - conf) / 2),
                                      ceiling((ndraw + 1) * (1 + conf) / 2))]))
}
