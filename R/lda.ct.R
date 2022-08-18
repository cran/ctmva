#'  Continuous-time Fisher's linear discriminant analysis
#'
#' A continuous-time version of Fisher's LDA, in which segments of the time
#' interval take the place of groups of observations.
#'
#' @param fdobj  continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param partition  a priori break points dividing the time interval into segments
#' @return  Object of class "\code{lda.ct}", a list consisting of
#' \item{scaling }{matrix of coefficients defining the discriminants (as in \code{\link[MASS]{lda}})}
#' \item{values }{eigenvalues giving the ratios of between to within sums of squares}
#' \item{partition }{the supplied \code{partition}}
#' \item{fdobj }{linear discriminants represented as an \code{"\link[fda]{fd}"} object}
#' \item{nld }{number of linear discriminants}
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso  \code{\link{plot.lda.ct}}; \code{\link[MASS]{lda}}, for the classical version
#' @examples
#'
#' ## see end of example in ?pca.ct
#'
#' @export lda.ct
lda.ct <-
  function(fdobj, partition) {
    bsb <- fdobj$basis
    C <- fdobj$coef
    allpts <- sort(c(bsb$rangeval, partition))
    nparts <- length(allpts) - 1
    nld <- min(length(partition), ncol(C))
    W. <- 0
    for (j in 1:nparts) {
      W. <- W. + inprod.cent(bsb,rng=c(allpts[j],allpts[j+1]))
    }
    svdW <- svd(t(C) %*% W. %*% C)
    disr <- svdW$d
    disr[disr<0] <- 0
    disr[disr>0] <- 1/sqrt(disr[disr>0])
    uv <- (svdW$u + svdW$v) / 2
    Wisr <- uv %*% diag(disr) %*% t(uv)
    B. <- inprod.cent(bsb) - W.
    B <- t(C) %*% B. %*% C
    W..BW.. <- Wisr %*% B %*% Wisr
    W..BW.. <- (W..BW.. + t(W..BW..)) / 2
    eig <- eigen(W..BW..)
    lst <- list(scaling = Wisr %*% eig$vectors[,1:nld],
                values = eig$values[1:nld],
                partition = partition,
                nld = nld,
                fdobj = fd(coef=C %*% eig$vectors[ , 1:nld], basisobj=bsb, 
                           fdnames=paste0("LD",1:nld)))
    class(lst) <- "lda.ct"
    lst
  }
