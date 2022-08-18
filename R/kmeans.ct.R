#'  Continuous-time k-means clustering
#'
#' A continuous-time version of k-means clustering in which each clusters is a
#' time segments or set of time segments.
#'
#' @param fdobj continuous-time multivariate data set of class \code{"\link[fda]{fd}"}
#' @param k  number of clusters
#' @param common_trend  logical: Should the curves be centered with respect to the mean function?
#'  Defaults to \code{FALSE}.
#' @param init.pts  a set of k time points. The observations at these time points
#' serve as initial values for the k means. Randomly generated if not supplied.
#' @param tol  convergence tolerance for the k means
#' @param max.iter  maximum number of iterations
#' @return  Object of class "\code{kmeans.ct}", a list consisting of
#' \item{fdobj }{the supplied \code{fdobj}}
#' \item{means }{means of the k clusters}
#' \item{transitions }{transition points between segments}
#' \item{cluster }{cluster memberships in the segments defined by the transitions}
#' @author Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#'
#' @seealso  \code{\link{plot.kmeans.ct}}
#' @examples
#'
#'
#' require(fda)
#' data(CanadianWeather)
#' daybasis <- create.bspline.basis(c(0,365), nbasis=55)
#' tempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis)$fd
#' kmtemp3 <- kmeans.ct(tempfd, 3)
#' plot(kmtemp3)
#'
#' @importFrom polynom poly.calc
#' @importFrom stats coef
#'
#' @export kmeans.ct
kmeans.ct <-
  function(fdobj, k, common_trend=FALSE, init.pts=NULL, tol=.001, max.iter=100) {
    if (common_trend) fdobj <- fda::center.fd(fdobj)

    # Select initial m_1,...,m_k
    bsb <- fdobj$basis
    if (bsb$type != "bspline") stop("Continuous-time k-means implemented only for B-spline bases.")
    rng <- bsb$rangeval
    if (is.null(init.pts)) init.pts <- sort(runif(k, rng[1], rng[2]))
    means <- eval.fd(init.pts, fdobj)	# k x p
    oldmeans <- means + 2*tol

    norder <- bsb$nbasis - length(bsb$params)  # = 4
    breaks <- c(rng[1], bsb$params, rng[2])
    nintervals <- length(breaks) - 1

    B <- t(fdobj$coef)    # p x nbasis
    count <- 1

    # Obtain coefficients of polynomial within interval
    oldoption <- options()  
    on.exit(options(oldoption))
    options(digits=15)
    scalefac <- 10^floor(log10(rng[2]))
    polycoefs <- array(dim=c(norder, bsb$nbasis, nintervals))
    for (int in 1:nintervals) {
      x <- seq(breaks[int], breaks[int+1], length = norder)
      xx <- x / scalefac    # rescale x to xx to avoid underflow in polynomial coefficients
      Phi <- eval.basis(x, bsb)
      polycoefs[ , , int] <- as.matrix(data.frame(apply(Phi, 2, function(v) coef(polynom::poly.calc(xx,v)))))
    }
    rawcoefs <- polycoefs
    for (nn in 2:norder) polycoefs[nn,,] <- polycoefs[nn,,] / scalefac^(nn-1)  # rescale back

    while (mean(abs(means-oldmeans)) > tol && count<=max.iter) {
      message("Iteration", count, "\n")
      oldmeans <- means

      # Find all crossing points of distance-from-cluster-center functions
      crosspts <- c()
      for (i in 1:(k-1)) for (j in (i+1):k) {
        for (int in 1:nintervals) {
          # Find zeroes of polynomial equal to d_i(t)^2 - d_j(t)^2 within subinterval
          coefvec <- 2 * polycoefs[,,int] %*% crossprod(B, means[j,]-means[i,]) + c(sum(means[i,]^2-means[j,]^2),0,0,0)
          zeros <- polyroot(coefvec)
          # If any of those zeros fall in the subinterval, add to crosspts
          realzeros <- Re(zeros[abs(Im(zeros)/Re(zeros))<1e-5])  # ad hoc criterion for a zero to be real...
          zeros_here <- realzeros[realzeros >= breaks[int] & realzeros<=breaks[int+1]]
          crosspts <- c(crosspts, zeros_here)
        }
      }
      crosspts <- sort(crosspts)

      # Determine minimum-distance cluster
      # for each "section" (= interval between consecutive crossing points)
      # --> Update cluster partition
      nsections <- length(crosspts) + 1
      testdists <- matrix(NA, nsections, k)
      starts <- c(rng[1], crosspts)
      ends <- c(crosspts, rng[2])
      testpts <- (starts + ends) / 2
      for (m in 1:k) testdists[ , m] <- apply(eval.fd(testpts, fdobj) - rep(1,nsections) %o% means[m,], 1, function(v) sum(v^2))
      cluster_memb <- apply(testdists, 1, which.min)

      # Update cluster means
      for (m in 1:k) {
        which_sections <- which(cluster_memb==m)
        cluster_length <- sum((ends - starts)[which_sections])

        # Integrate f over cluster-m sections
        int_f <- rep(0,ncol(fdobj$coef))
        for (section in which_sections) {

          # Integrate f within section (check it!)
          start. <- starts[section]; end. <- ends[section]
          pts <- c(start., breaks[breaks>start. & breaks<end.], end.)
          npts <- length(pts)
          midpts <- cbind(pts[-1],pts[-npts]) %*% rbind(1:5, 5:1) / 6
          allpts0 <- c(pts, as.vector(midpts))
          if (npts>2) allpts0 <- c(allpts0, pts[2:(npts-1)])
          allpts <- sort(allpts0)
          mult <- rep(diff(pts), each=7) * rep(c(41/140, 54/35, 27/140, 68/35, 27/140, 54/35, 41/140),(npts-1)) / 6
          int_section <- colSums(diag(mult) %*% eval.fd(allpts, fdobj))
          int_f <- int_f + int_section
        }
        means[m, ] <- int_f / cluster_length
      }

      count <- count+1
    }

    transitions <- ends[-length(ends)][diff(cluster_memb)!=0]
    cluster <- cluster_memb[c(1, 1+which(diff(cluster_memb)!=0))]
    lst <- list(fdobj=fdobj, means=means, transitions=transitions, cluster=cluster)
    class(lst) <- "kmeans.ct"
    lst
  }
