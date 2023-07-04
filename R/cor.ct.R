#' Continuous-time correlation or cross-correlation matrix
#'
#' Computes the correlation matrix of a continuous-time multivariate
#' data set represented as an \code{\link[fda]{fd}} object; or the cross-correlation
#' matrix of two such data sets.
#'
#' @param fdobj1  continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param fdobj2  an optional second data set
#' @param common_trend  logical: centering wrt mean function if \code{TRUE},
#' without centering if \code{FALSE} (the default)
#' @return  A matrix of (cross-) correlations
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso
#' \code{\link[fda]{center.fd}}, for centering of \code{"\link[fda]{fd}"} objects;  \code{\link{inprod.cent}}
#' @examples
#'
#'
#' # Canadian temperature data
#'
#' require(fda)
#' require(corrplot)
#' data(CanadianWeather)
#' daybasis <- create.fourier.basis(c(0,365), nbasis=55)
#' tempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis)$fd
#'
#' ## The following yields a matrix of correlations that are all near 1:
#' rawcor <- cor.ct(tempfd)
#' corrplot(rawcor, method = 'square', type = 'lower', tl.col="black", tl.cex = 0.6)
#' ## This occurs due to a strong seasonal trend that is common to all stations
#' ## Removing this common trend leads to a more interesting result:
#' dtcor  <- cor.ct(tempfd, common_trend = TRUE)
#' ord <- corrMatOrder(dtcor)
#' dtcord <- dtcor[ord,ord]
#' corrplot(dtcord, method = 'square', type = 'lower', tl.col="black", tl.cex = 0.6)
#'
#'
#' @export cor.ct
cor.ct <-
  function(fdobj1, fdobj2=fdobj1, common_trend=FALSE) {
    if(common_trend){
      fdobj1 <- fda::center.fd(fdobj1)
      fdobj2 <- fda::center.fd(fdobj2)
    }
    P0 <- inprod.cent(fdobj1$basis, fdobj2$basis)
    Phi1 <- t(fdobj1$coef)
    Phi2 <- t(fdobj2$coef)
    covmat <- Phi1 %*% P0 %*% t(Phi2) / diff(fdobj1$basis$rangeval)
    dinv1 <- matrix(0,nrow(Phi1),nrow(Phi1))
    diag(dinv1) <- 1/sqrt(diag(Phi1 %*% P0 %*% t(Phi1) / diff(fdobj1$basis$rangeval)))
    dinv2 <- matrix(0,nrow(Phi2),nrow(Phi2))
    diag(dinv2) <- 1/sqrt(diag(Phi2 %*% P0 %*% t(Phi2) / diff(fdobj2$basis$rangeval)))
    matt <- dinv1 %*% covmat %*% dinv2
    matt[matt > 1] <- 1
    matt[matt < -1] <- -1

    if(!is.null(fdobj1$fdnames$reps)) rownames(matt) <- paste0(fdobj1$fdnames$reps)
    else rownames(matt) <- paste0("Comp.", 1:ncol(fdobj1$coef))

    if(!is.null(fdobj2$fdnames$reps)) colnames(matt) <- paste0(fdobj2$fdnames$reps)
    else colnames(matt) <- paste0("Comp.", 1:ncol(fdobj2$coef))
    return(matt)
}
