#'  Continuous-time canonical correlation analysis
#'
#'  A continuous-time version of canonical correlation analysis (CCA).
#'
#' @param fdobj1,fdobj2  a pair of continuous-time multivariate data sets, of class \code{"\link[fda]{fd}"}
#' @return  A list consisting of
#' \item{vex1, vex2 }{matrices defining the canonical variates. The first columns of each give the coefficients defining the first pair of canonical variates; and so on.}
#' \item{cor }{canonical correlations, i.e., correlations between the pairs of canonical variates}
#'
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @note  Columns of the output matrix \code{vex2} are flipped as needed to ensure positive correlations.
#' @seealso  \code{\link{cancor}}, for classical CCA
#' @examples
#'
#'
#' # CCA relating Canadian daily temperature and precipitation data
#' require(fda)
#' data(CanadianWeather)
#' daybasis <- create.bspline.basis(c(0,365), nbasis=80)
#' tempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis)$fd
#' precfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"log10precip"], daybasis)$fd
#' tpcor <- cca.ct(tempfd, precfd)
#' par(mfrow=1:2)
#' barplot(tpcor$vex1[,1], horiz=TRUE, las=1, main="Temperature",
#'             sub="First canonical coefficients vector")
#' barplot(tpcor$vex2[,1], horiz=TRUE, las=1, main="Log precipitation",
#'             sub="First canonical coefficients vector")
#'
#'
#'
#' @export cca.ct
cca.ct <-
  function(fdobj1, fdobj2) {
    if (!isTRUE(fdobj1$basis == fdobj2$basis)) stop("Basis for both fd arguments must be same.")
    if (fdobj1$basis$type != "bspline")
      stop("Continuous-time canonical correlation is implemented only for B-spline bases.")
    if (fdobj1$basis$nbasis < ncol(fdobj1$coefs) + ncol(fdobj2$coefs))
      stop("Combined number of functions in fdobj1, fdobj2 cannot exceed the number of basis functions.")
    P0 <- inprod.cent(fdobj1$basis)

    C1 <- fdobj1$coef
    C2 <- fdobj2$coef
    if (qr(crossprod(C1, P0%*%C1))$rank!=ncol(fdobj1$coefs))
      stop("Covariance matrix for fdobj1 is not full-rank.")
    i11 <- solve(crossprod(C1, P0%*%C1))
    m12 <- crossprod(C1, P0%*%C2)
    if (qr(crossprod(C2, P0%*%C2))$rank!=ncol(fdobj2$coefs))
      stop("Covariance matrix for fdobj2 is not full-rank.")
    i22 <- solve(crossprod(C2, P0%*%C2))
    mat1 <- i11 %*% m12 %*% i22 %*% t(m12)
    mat2 <- i22 %*% t(m12) %*% i11 %*% m12
    eig1 <- eigen(mat1); eig2 <- eigen(mat2)
    vex1 <- eig1$vec
    vex2 <- eig2$vec
    #if (any(is.complex(vex1[,min(ncol(vex1),ncol(vex2))]))||any(is.complex(vex2[,min(ncol(vex1),ncol(vex2))])))
    #  warning("Complex canonical coefficients. These may correspond to perfectly correlated canonical variates.")
    rownames(vex1) <- colnames(fdobj1$coef)
    rownames(vex2) <- colnames(fdobj2$coef)

    # Flip signs to get positive correlations
    for (i in 1:min(ncol(vex1),ncol(vex2))) {
      if (crossprod(fdobj1$coef%*%Re(vex1[,i]), P0 %*% fdobj2$coef%*%Re(vex2[,i])) < 0)
        vex2[,i] <- -vex2[,i]
    }
    list(vex1=vex1, vex2=vex2, cor=sqrt((eig1$val)[1:min(ncol(vex1),ncol(vex2))]))
  }
