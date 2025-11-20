#' Credible interval estimation for curve correlation based on posterior simulation
#'
#' Inputs raw data representing two curves, and computes a credible interval
#' for the curve correlation between them simulating from the approximate posterior
#' distribution of the joint spline coefficient vector.
#'
#' @param y,time,curve,k,min.overlap see \code{\link{ccor}}
#' @param method \code{"indep"} (curves fitted independently)
#' or \code{"mvn"} (curves fitted together, with error correlation estimated based on multivariate normal likelihood)
#' @param conf confidence level
#' @param ndraw number of draws from posterior distribution of spline coefficient vector
#' @return A list with components
#' \item{cor }{curve correlation}
#' \item{model }{the model for the two curves (if \code{method=="mvn"}), or a list of the two curves' models  (if \code{method=="indep"})}
#' \item{bsb }{B-spline basis (from package \code{fda}) used for the curves}
#' \item{Vc.fda }{corrected posterior covariance matrix for the coefficients with respect to the B-spline basis \code{bsb}
#' (see the component \code{$Vc} in \code{\link[mgcv]{gamObject}})}
#' \item{sims }{curve correlations for the \code{ndraw} draws from the posterior distribution}
#' \item{ci }{credible interval for the curve correlation}
#'
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>, Noemi Foa, Dror Arbiv and Biplab Paul <paul.biplab497@gmail.com>
#'
#'
#' @seealso \code{\link{ccor}}, \code{\link{ccor_boot}}, \code{\link[mgcv]{mvn}}
#' @examples
#'
#' ## Not run:
#' if (interactive () &&
#'     requireNamespace("wbwdi", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE)) {
#'
#'     # Curve correlation of per capita GDP and fertility rate in Paraguay
#'     wdi_dat <- wbwdi::wdi_get(entities = c("PRY"), start_year=1960, end_year=2023,
#'                        indicators = c("NY.GDP.PCAP.KD","SP.DYN.TFRT.IN"), format="wide") |>
#'         dplyr::rename(percapitaGDP = NY.GDP.PCAP.KD, fertility = SP.DYN.TFRT.IN)
#'     ggplot2::ggplot(wdi_dat, aes(percapitaGDP, fertility, color=year)) + geom_point()
#'
#'     y <- as.matrix(wdi_dat[ , c("percapitaGDP", "fertility")])
#'
#'     set.seed(345)
#'
#'     ci <- list()
#'     ci[[1]] <- ccor_posim(y=y, time=wdi_dat$year, method="indep")
#'     ci[[2]] <- ccor_posim(y=y, time=wdi_dat$year, method="mvn")
#'     ci[[3]] <- ccor_boot(y=y, time=wdi_dat$year, ndraw=399)
#'
#'     tabl <- matrix(NA, 3, 3)
#'     for (k in 1:3) tabl[k, ] <- c(ci[[k]]$cor, ci[[k]]$ci)
#'     dimnames(tabl) <- list(c("Posim_indep", "Posim_MVN", "Bootstrap"), c("Est","Lower95","Upper95"))
#'     round(tabl, 4)
#' }
#' ## End(Not run)
#'
#' @importFrom MASS ginv
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom mgcv gam
#' @importFrom mgcv mvn
#' @importFrom fda fd
#' @importFrom fda create.bspline.basis
#' @importFrom stats fitted
#' @importFrom stats model.matrix
#' @export

ccor_posim <-
function(y, time, curve= NULL, method, k=15, conf=.95, ndraw=999, min.overlap=0) {

    # Check that time is a vector
    if (!is.vector(time)) stop("time must be a vector")

    # Check that y is either a vector or a 2-column matrix
    if (!is.vector(y) && !(is.matrix(y) && ncol(y) == 2)) {
        stop("y must be either a vector or a 2-column matrix")
    }

    if (is.vector(y)) {
        if (method=="mvn") stop("If y is a vector, method must be 'indep'")
        if (length(time) != length(y) || length(curve) != length(y)) {
            stop("If y is a vector, time and curve must be vectors of the same length as y")
        }
        curve <- factor(curve)
        levs <- levels(curve)
        if (length(levs) != 2)
            stop("If y is a vector, curve must take exactly two values")
        ylist <- list(y[curve == levs[1]], y[curve == levs[2]])
    }

    if (is.matrix(y)) {
        if (nrow(y) != length(time))
            stop("If y is a matrix, its number of rows must match the length of time")
        ylist <- list(y[, 1], y[, 2])
    }

    if (method == 'indep') {
        model <- wots <- nots <- bsb <- ginvB <- nucoef <-
            est_y <- Vc.mgcv <- TT <- Vc.fda <- simC <- list()
        for (j in 1:2) {
            t_j <- if (is.vector(y)) time[curve == levs[j]] else time
            model[[j]] <- gam(ylist[[j]] ~ s(t_j, bs='bs', k=k), method="REML")
            wots[[j]] <- model[[j]]$smooth[[1]]$knots
            nots[[j]] <- sort(c(range(t_j), wots[[j]][wots[[j]]>min(t_j) & wots[[j]]<max(t_j)]))
            bsb[[j]] <- create.bspline.basis(range(t_j), breaks=nots[[j]])
            ginvB[[j]] <- ginv(fda::eval.basis(t_j, bsb[[j]]))
            nucoef[[j]] <- ginvB[[j]] %*% fitted(model[[j]])
            est_y[[j]] <- fd(coef = nucoef[[j]], basisobj = bsb[[j]])
            Vc.mgcv[[j]] <- model[[j]]$Vc
            TT[[j]] <- ginvB[[j]] %*% model.matrix(model[[j]])
            Vc.fda[[j]] <- nearPD(TT[[j]] %*% Vc.mgcv[[j]] %*% t(TT[[j]]))$mat
            simC[[j]] <- mvrnorm(ndraw, mu=as.vector(nucoef[[j]]), Sigma=Vc.fda[[j]])
        }
        overlap <- c(max(bsb[[1]]$rangeval[1], bsb[[2]]$rangeval[1]),
                     min(bsb[[1]]$rangeval[2], bsb[[2]]$rangeval[2]))
        if (diff(overlap) < min.overlap) stop("Insufficient overlap of time ranges")
        Q12 <- inprod.cent(bsb[[1]], bsb[[2]]) / diff(overlap)
        Q1 <- inprod.cent(bsb[[1]], rng = overlap) / diff(overlap)
        Q2 <- inprod.cent(bsb[[2]], rng = overlap) / diff(overlap)
        v12 <- rowSums((simC[[1]] %*% Q12) * simC[[2]])
        v1  <- rowSums((simC[[1]] %*% Q1) * simC[[1]])
        v2 <-  rowSums((simC[[2]] %*% Q2) * simC[[2]])
    }
    else if (method == 'mvn') {
        y1 <- ylist[[1]]; y2 <- ylist[[2]]
        model <- gam(list(y1~s(time, bs='bs', k=k),
                          y2~s(time, bs='bs', k=k)), family=mvn(2))
        wots <- model$smooth[[1]]$knots
        nots <- sort(c(range(time), wots[wots>min(time) & wots<max(time)]))
        bsb <- create.bspline.basis(range(time), breaks=nots)
        ginvB <- ginv(fda::eval.basis(time, bsb))
        nucoef <- ginvB %*% fitted(model)
        est_y1 <- fd(coef = nucoef[,1], basisobj = bsb)
        est_y2 <- fd(coef = nucoef[,2], basisobj = bsb)
        Vc.mgcv <- model$Vc[1:(2*k), 1:(2*k)]
        TT <- matrix(0, 2*k, 2*k)
        TT[1:k, 1:k] <- ginvB %*% model.matrix(model)[,1:k]
        TT[(k+1):(2*k), (k+1):(2*k)] <- ginvB %*% model.matrix(model)[,(k+1):(2*k)]
        Vc.fda <- nearPD(TT %*% Vc.mgcv %*% t(TT))$mat
        simC <- mvrnorm(ndraw, mu=as.vector(nucoef), Sigma=Vc.fda)
        simC1 <- simC[, 1:k]
        simC2 <- simC[, (k+1):(2*k)]
        Q <- inprod.cent(bsb)/diff(range(time))
        simC1Q <- simC1 %*% Q
        v12 <- rowSums(simC1Q * simC2)
        v1  <- rowSums(simC1Q * simC1)
        v2 <-  rowSums((simC2 %*% Q) * simC2)
    }

    ctcor.vec <- v12/(sqrt(v1*v2))

    list(cor = ifelse(method == 'mvn',
                      cor.ct(est_y1, est_y2), cor.ct(est_y[[1]], est_y[[2]])),
         model=model,
         bsb=bsb,
         Vc.fda=Vc.fda,
         sims=ctcor.vec,
         ci=sort(ctcor.vec)[c(floor((ndraw+1)*(1-conf)/2),ceiling((ndraw+1)*(1+conf)/2))]
    )
}
