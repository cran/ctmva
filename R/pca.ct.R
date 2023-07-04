#'  Continuous-time principal component analysis
#'
#'  A continuous-time version of principal component analysis.
#'
#' @param fdobj continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param cor  logical: use correlation matrix if \code{TRUE},
#' covariance if \code{FALSE} (the default)
#' @param common_trend  logical: Should the curves be centered with respect to the mean function?
#'  Defaults to \code{FALSE}.
#' @return Returns a list including:
#' \item{var}{variances of the principal components.}
#' \item{loadings}{the matrix of loadings (i.e., its columns are the
#' eigenvectors of the continuous-time covariance).}
#'  \item{scorefd}{score functions.}
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso  \code{\link{cov.ct}}; \code{\link{princomp}}, for the classical version
#' @examples
#'
#'
#' # Data for one session from a classic EEG data set
#' require(fda)
#' require(eegkit)
#' data(eegdata)
#' data(eegcoord)
#' longdat <- subset(eegdata, subject=="co2a0000369" & trial==0)
#' widedat <- reshape(longdat, direction="wide", drop=c("subject","group","condition","trial"),
#'                  v.names="voltage",idvar="channel")
#'
#' # Convert time series for 64 channels to a functional data object
#' bsb <- create.bspline.basis(c(0,255),nbasis=30)
#' fdo <- Data2fd(argvals=0:255, y=t(as.matrix(widedat[,-1])), basisobj=bsb)
#' plot(fdo)
#'
#' # Now do PCA and display first loadings for 3 PC's,
#' # along with percent variance explained by each
#' pcc <- pca.ct(fdo)
#' pve <- 100*pcc$var/sum(pcc$var)
#' par(mfrow=c(1,3))
#' cidx <- match(widedat[,1],rownames(eegcoord))
#' eegspace(eegcoord[cidx,4:5],pcc$loadings[,1], colorlab="PC1 loadings",
#'          main=paste0(round(pve[1],0), "%"), mar=c(17,3,12,2), cex.main=2)
#' eegspace(eegcoord[cidx,4:5],pcc$loadings[,2], colorlab="PC2 loadings",
#'          main=paste0(round(pve[2],0), "%"), mar=c(17,3,12,2), cex.main=2)
#' eegspace(eegcoord[cidx,4:5],pcc$loadings[,3], colorlab="PC3 loadings",
#'          main=paste0(round(pve[3],0), "%"), mar=c(17,3,12,2), cex.main=2)
#'
#' # Linear discriminant analysis: discriminating among the 1st, 2nd and 3rd portions
#'#  of the time interval
#' ld <- lda.ct(fdo, c(85,170))
#' plot(ld)
#' eegspace(eegcoord[cidx,4:5],ld$scaling[,1], colorlab="LD1 coefficients",
#'          mar=c(17,3,12,2), cex.main=2)
#' eegspace(eegcoord[cidx,4:5],ld$scaling[,2], colorlab="LD2 coefficients",
#'          mar=c(17,3,12,2), cex.main=2)
#'
#'
#' @export pca.ct
pca.ct <- function(fdobj, cor=FALSE, common_trend=FALSE) {
    if (common_trend) fdobj <- fda::center.fd(fdobj)
    covmat <- cov.ct(fdobj)
    if (cor) {
        isr <- diag(1/sqrt(diag(covmat)))
        cmat <- isr %*% covmat %*% isr
    } else cmat <- covmat
    output <- eigen(cmat)
    colnames(output$vectors) <- paste0("Comp.", 1:ncol(fdobj$coef))
    rownames(output$vectors) <- fdobj$fdnames[[2]]
    results <- list(var = output$values, loadings = output$vectors)
    ccoef <- center.ct(fdobj)$coefs
    if (cor) ccoef <- ccoef %*% isr
    results$scorefd <- fd(ccoef %*% output$vectors, fdobj$basis)
    return(results)
}
