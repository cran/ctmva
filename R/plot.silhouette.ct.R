#' Plot a silhouette.ct object
#'
#' Plots the silhouette index, generated by a call to
#' \code{\link{silhouette.ct}}, for a continuous-time k-means clustering object.
#' 
#' @param x silhouette object produced by \code{\link{silhouette.ct}}
#' @param mark.transitions logical: Should transitions between clusters be marked
#' with vertical lines? Defaults to \code{TRUE}.
#' @param xlab,ylab  x- and y-axis labels
#' @param \dots other arguments passed to \code{\link{plot}}
#' @return None; a plot is generated.
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#' @seealso \code{\link{kmeans.ct}}, which includes an example; \code{\link{silhouette.ct}}
#' @export plot.silhouette.ct
#' @export
plot.silhouette.ct <- function(x, mark.transitions=TRUE, xlab="Time", ylab="Silhouette", ...) {
    with(x, plot(grid, value, type='l', xlab=xlab, ylab=ylab, ...))
    if (mark.transitions) {
        for (pt in x$transitions) abline(v=pt, col="grey")
        midpts <- (c(min(x$grid),x$transitions) + c(x$transitions,max(x$grid))) / 2
        mtext(x$cluster, at=midpts)
    }
}
