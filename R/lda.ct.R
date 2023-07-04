#'  Continuous-time Fisher's linear discriminant analysis
#'
#' A continuous-time version of Fisher's LDA, in which segments of the time
#' interval take the place of groups of observations.
#'
#' The \code{means} and \code{scaling} components of the output are similar to
#' \code{\link[MASS]{lda}}, but unlike that function, \code{lda.ct} performs only
#' \emph{Fisher's} LDA and cannot incorporate priors or perform classification.
#'
#' @param fdobj  continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param partition  a priori break points dividing the time interval into segments
#' @param part.names  optional character vector of names for the segments
#' @return  Object of class "\code{lda.ct}", a list consisting of
#' \item{means }{means of the variables within each segment}
#' \item{scaling }{matrix of coefficients defining the discriminants (as in \code{\link[MASS]{lda}})}
#' \item{values }{eigenvalues giving the ratios of between to within sums of squares}
#' \item{partition }{the supplied \code{partition}}
#' \item{fdobj }{linear discriminants represented as an \code{"\link[fda]{fd}"} object}
#' \item{nld }{number of linear discriminants}
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso  \code{\link{plot.lda.ct}}; \code{\link[MASS]{lda}}, for the classical version
#' @examples
#'
#' ## see end of example in ?pca.ct
#'
#' @export lda.ct
lda.ct <-
  function(fdobj, partition, part.names=NULL) {
    bss <- fdobj$basis
    C <- fdobj$coef
    allpts <- sort(c(bss$rangeval, partition))
    nparts <- length(allpts) - 1
    nld <- min(length(partition), ncol(C))
    W. <- 0
    meanbasis.parts <- matrix(0,nparts,bss$nbasis)
    for (j in 1:nparts) {
      W. <- W. + inprod.cent(bss,rng=c(allpts[j],allpts[j+1]))
      meanbasis.parts[j,] <- meanbasis(bss, rng=c(allpts[j],allpts[j+1]))
    }
    svdW <- svd(t(C) %*% W. %*% C)
    disr <- svdW$d
    disr[disr<0] <- 0
    disr[disr>0] <- 1/sqrt(disr[disr>0])
    uv <- (svdW$u + svdW$v) / 2
    Wisr <- uv %*% diag(disr) %*% t(uv)
    B. <- inprod.cent(bss) - W.
    B <- t(C) %*% B. %*% C
    W..BW.. <- Wisr %*% B %*% Wisr
    W..BW.. <- (W..BW.. + t(W..BW..)) / 2
    eig <- eigen(W..BW..)
    scaling <- Wisr %*% eig$vectors[,1:nld]
    colnames(scaling) <- paste0("LD", 1:nld)
    means <- meanbasis.parts %*% C
    if (is.null(part.names)) part.names <- paste0("part",1:nparts)
    rownames(means) <- part.names
    rownames(scaling) <- colnames(means) <-  fdobj$fdnames[[2]]

    lst <- list(means = means,
                scaling = scaling,
                values = eig$values[1:nld],
                partition = partition,
                nld = nld,
                fdobj = fd(coef=C %*% scaling, basisobj=bss,
                                      fdnames=paste0("LD",1:nld)))
    class(lst) <- "lda.ct"
    lst
}
