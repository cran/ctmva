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
#' @return  If \code{fdobj1} and \code{fdobj2} each consist of a single curve,
#' the (scalar) CT correlation between them. Otherwise a matrix of (cross-) correlations is returned.
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @note When \code{fdobj1==fdobj2}, \code{cor.ct(\dots)} is the same as \code{cov2cor(cov.ct(\dots))}.
#' @seealso
#' \code{\link[fda]{center.fd}}, for centering of \code{"\link[fda]{fd}"} objects;  \code{\link{inprod.cent}}
#' @examples
#'
#' \dontrun{
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
#'}
#'
#' @export cor.ct
cor.ct <-
  function(fdobj1, fdobj2=fdobj1, common_trend=FALSE) {
    if (common_trend) {
      fdobj1 <- fda::center.fd(fdobj1)
      fdobj2 <- fda::center.fd(fdobj2)
    }
    rng1 <- fdobj1$basis$rangeval; rng2 <- fdobj2$basis$rangeval
    overlap <- c(max(rng1[1],rng2[1]), min(rng1[2],rng2[2]))
    P0 <- inprod.cent(fdobj1$basis, fdobj2$basis) # range of integration automatically set to 'overlap'
    Phi1 <- t(fdobj1$coef)
    Phi2 <- t(fdobj2$coef)
    covmat <- Phi1 %*% P0 %*% t(Phi2) / diff(overlap)
    dinv1 <- matrix(0,nrow(Phi1),nrow(Phi1))
    samebasis <- (length(all.equal(fdobj1$basis, fdobj2$basis))==1)
    rng. <- if (samebasis) NULL else overlap
    diag(dinv1) <- 1/sqrt(diag(Phi1 %*% inprod.cent(fdobj1$basis,rng=rng.) %*% t(Phi1) / diff(overlap)))
    dinv2 <- matrix(0,nrow(Phi2),nrow(Phi2))
    diag(dinv2) <- 1/sqrt(diag(Phi2 %*% inprod.cent(fdobj2$basis,rng=rng.) %*% t(Phi2) / diff(overlap)))
    matt <- dinv1 %*% covmat %*% dinv2
    matt[matt > 1] <- 1
    matt[matt < -1] <- -1

    if(!is.null(fdobj1$fdnames$reps)) rownames(matt) <- fdobj1$fdnames$reps
    else rownames(matt) <- paste0("Comp.", 1:ncol(fdobj1$coef))

    if(!is.null(fdobj2$fdnames$reps)) colnames(matt) <- fdobj2$fdnames$reps
    else colnames(matt) <- paste0("Comp.", 1:ncol(fdobj2$coef))

    if (length(matt)==1) matt <- matt[1,1]

    return(matt)
}
