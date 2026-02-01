#' Continuous-time multidimensional scaling
#'
#' A continuous-time version of classical multidimensional scaling.
#'
#' @param d matrix of dissimilarities
#' @param argvals time points
#' @param weights quadrature weights for the time points: either "equal" (default) 
#' or "trap" (trapezoidal, recommended for unequally spaced times)
#' @param nbasis number of B-spline basis functions
#' @param norder order of B-splines (degree plus 1); default is 4 (cubic splines)
#' @param breaks knots, including the endpoints of the time range
#' @param k number of principal coordinates
#' @param smooth bivariate smoothing method: "sandwich" (default) for the sandwich smoother,
#' "one-step" for a more traditional but slower tensor product smoother
#' @param lambda smoothing parameter
#' @param gcvplot logical: Should a plot of log lambda versus the GCV criterion be produced?
#' @param Im.tol tolerance for imaginary component of eigenvalues: imaginary components
#'  with magnitude below this is set to zero, those above it generate an error.
#' @param recenter logical: Should the solution be double-centered?
#'
#' @return An object of class "mds.ct", which is a list with components
#' \item{funcs}{an object of class "\code{\link[fda]{fd}}" giving the \eqn{k} principal
#' coordinate functions}
#' \item{evals}{eigenvalues}
#' \item{argvals}{the given time points}
#' \item{GOF}{a \eqn{k \times 2} matrix giving the proportion of dissimilarity explained, 
#' according to the two definitions used by the \code{GOF} component of the output from \code{\link{cmdscale}}}
#' \item{sandwich}{output of \code{\link{s3}}, if \code{smooth=="sandwich"}}
#' \item{recenter}{value of the argument \code{recenter}}
#' 
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il> 
#'
#' @seealso \code{\link{cmdscale}}, for the classical version
#' @importFrom fda create.bspline.basis eval.basis fd getbasispenalty
#' @importFrom mgcv gam
#' @export
#'
#' @examples
#'
#' if (interactive() && 
#'     requireNamespace("fda", quietly = TRUE) &&
#'     requireNamespace("viridisLite", quietly = TRUE) &&
#'     requireNamespace("vegan", quietly = TRUE)) {
#'     
#'     require(ggplot2)
#'     
#'     data("handwrit", package = "fda")
#'     
#'     fives <- 5 * 0:280  + 1
#'     hw <- handwrit[ , 1, ]
#'     sd. <- .015
#'     
#'     hh <- cbind(hw[fives,], rnorm(281, sd=sd.))
#'     classical <- cmdscale(dist(hh), eig=TRUE)
#'     ctmds <- mds.ct(dist(hh), argvals=handwritTime[fives], nbasis=100, recenter=TRUE,
#'                     weights="trap", norder=7, lambda=exp(2:40))  
#'                     # a plot of GCV versus lambda is produced
#'     
#'     pro.classical <- vegan::procrustes(hw[fives,], classical, scale=FALSE)
#'     dat.classical <- as.data.frame(pro.classical$Yrot)
#'     pro.ctmds <- vegan::procrustes(hw[fives,], fda::eval.fd(ctmds$argvals, ctmds$funcs),
#'                                    scale=FALSE)
#'     dat.ctmds <- as.data.frame(matrix(pro.ctmds$translation, length(handwritTime), 2, byrow=TRUE) +
#'                                    fda::eval.fd(handwritTime, ctmds$funcs) %*% pro.ctmds$rotation)
#'     names(dat.classical) <- names(dat.ctmds) <- c("x", "y")
#'     dat.classical$time <- handwritTime[fives]
#'     dat.ctmds$time <- handwritTime
#'     
#'     # Plot of classical (dots) versus continuous-time (curve) MDS reconstructions
#'     # of the handwritten "fda" (grey), corrupted by noise
#'     g1 <- ggplot(hw, aes(x=X, y=Y)) + geom_point(color="darkgrey", size=.6) + coord_fixed() +
#'         geom_point(data = dat.classical, aes(x, y, color = time), size=1) +
#'         geom_point(data = dat.ctmds, aes(x, y, color = time), size=.6) +
#'         scale_color_gradientn(colors=viridisLite::plasma(50)) + 
#'         labs(x="", y="", color="Time (ms)", 
#'              title=bquote(sigma == .(sd.) * ", " 
#'                           ~ .(round(100*ctmds$GOF[2,1], 1)) * "% explained")) 
#'     print(g1)
#'     
#'     g2 <- plot(ctmds)  # uses plot.mds.ct
#'     print(g2)
#'     
#'     g3 <- procrustes.plot.mds(ctmds, dat.classical[ , 1:2])
#'     print(g3)
#'     
#' }
#' 
mds.ct <-
    function(d, argvals=NULL, weights="equal", nbasis=10, norder=4, breaks=NULL, k=2, 
             smooth="sandwich", lambda=exp(-10:20), gcvplot=TRUE, Im.tol=1e-15, recenter=TRUE) {
        if (is.null(lambda)) stop("lambda must be given")
        if (is.null(argvals)) argvals <- as.numeric(rownames(d))
        basis <- create.bspline.basis(range(argvals), nbasis, norder, breaks)
        if (!is.matrix(d)) d <- as.matrix(d)
        if (anyNA(d) || any(!is.finite(d))) stop("d contains NA or non-finite values")
        if (!isSymmetric(d, tol = 1e-12)) d <- (d + t(d)) / 2
        diag(d) <- 0
        n <- nrow(d)
        if (diff(range(diff(argvals)))>0 && weights=="equal") 
            warning("Unequally-spaced argvals; trapezoidal weights might be more appropriate") 
        if (weights == "equal") wtvec <- rep(1/n,n)
        else if (weights == "trap") wtvec <- trapwts(argvals)
        A <- -d^2/2
        v_row <- as.numeric(wtvec %*% A)    # 1 x n 
        term1 <- matrix(rep(v_row, n), nrow = n, byrow = TRUE)  # n x n 
        v_col <- as.numeric(A %*% wtvec)    # n x 1 
        term2 <- matrix(rep(v_col, n), nrow = n, byrow = FALSE)  # n x n 
        s_val <- as.numeric(wtvec %*% (A %*% wtvec))  # scalar
        Bmat <- A - term1 - term2 + s_val
        Bmat <- (Bmat + t(Bmat)) / 2
        M <- getbasispenalty(basis,0)
        P <- getbasispenalty(basis, norder-2)
        B <- eval.basis(argvals, basis)
        if (smooth=="one-step") {
            evalmat <- B %x% B
            Bmod <- gam(as.vector(Bmat) ~ evalmat - 1, paraPen=list(evalmat=list(M%x%P + P%x%M)), method="REML")
            coefmat <- matrix(Bmod$coef, basis$nbasis)
        } else if (smooth=="sandwich") {
            wich <- s3(Y=Bmat, B=B, P=P, lambda=lambda)
            coefmat <- wich$coefmat
            if (gcvplot) with(wich, 
                              plot(as.numeric(names(gcv)), gcv, log="x",
                                   xlab=expression(lambda), ylab="GCV"))
        }
        R <- (coefmat + t(coefmat)) / 2
        if (recenter) {
            mb <- meanbasis(basis)
            Rmb <- as.numeric(R %*% mb)
            Rmbmat <- matrix(rep(Rmb,nbasis), nbasis)
            R <- R - Rmbmat - t(Rmbmat) + matrix(as.numeric(crossprod(mb, Rmb)), nbasis, nbasis)
        }
        eig <- mygen(R %*% M, Im.tol=Im.tol)    
        evals <- eig$val
        if (max(abs(Im(evals))) > 1e-15) stop("Complex eigenvalues") else evals <- Re(evals)
        if (evals[k] < 0) stop(paste("First", k, "eigenvalues are not all positive")) 
        V <- eig$vec[,1:k]
        if (max(abs(Im(V))) > 1e-15) stop("Complex eigenvectors") else V <- Re(V)
        squared_norms_of_efuncs <- colSums(V * (M %*% V))
        GOF <- matrix(cumsum(evals[1:k]), k, 2)
        GOF[,1] <- GOF[,1] / sum(abs(evals))
        GOF[,2] <- GOF[,2] / sum(pmax(evals,0))
        dimnames(GOF) <- list(1:k, c("GOF(abs)", "GOF(pos.part)"))
        mdfd <- fd(coef = V %*% diag(sqrt(evals[1:k] / squared_norms_of_efuncs)), basisobj = basis)
        rslt <- list(funcs = mdfd, evals=evals, argvals=argvals, GOF=GOF,
                     sandwich=if (smooth=="sandwich") wich else NULL, recenter=recenter)
        class(rslt) <- "mds.ct"
        rslt
    }


#' Symmetric sandwich smoother 
#' 
#' This function is used by \code{\link{mds.ct}} to fit a symmetric version of 
#' the sandwich smoother of Xiao et al. (2013).
#'
#' @param Y a square symmetric matrix
#' @param B a B-spline basis evaluation matrix
#' @param P a penalty matrix
#' @param lambda a vector of smoothing parameter values
#' 
#' @return A list with components
#' \item{gcv}{generalized cross-validation (GCV) criterion for each value of \code{lambda}}
#' \item{bestlam}{value of \code{lambda} minimizing the GCV}
#' \item{coefmat}{matrix of tensor product B-spline coefficients}
#' \item{Yhat}{the smoothed version of \code{Y}}
#' 
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il> 
#'
#' @seealso \code{\link{mds.ct}} 
#' 
#' @references
#' Xiao, L., Li, Y. & Ruppert, D. (2013). Fast bivariate P-splines: the sandwich smoother.
#' \emph{Journal of the Royal Statistical Society Series B}, 75(3), 577–599.
#'
#' @export
#'
#' @examples
#' # see example for mds.ct (which calls s3)
s3 <- function(Y, B, P, lambda=exp(-10:20)) {
    if (is.null(lambda)) stop("lambda must be given")
    if (nrow(Y) != ncol(Y)) stop("Y must be square symmetric")
    if (max(abs(Y-t(Y))) > 1e-12)  {
        print(range(Y-t(Y)))
        stop("Y must be symmetric")
    }
    n <- nrow(Y)
    if (nrow(B) != n) stop("B must have same number of rows as Y")
    nlam <- length(lambda)
    BTB <- crossprod(B)
    svdBTB <- svd(BTB)
    if (any(svdBTB$d <= 0)) stop("B'B has negative singular values")
    BTBinvsqrt <- svdBTB$u %*% diag(1/sqrt(svdBTB$d)) %*% t(svdBTB$u)
    svd.. <- svd(BTBinvsqrt %*% P %*% BTBinvsqrt)
    A <- B %*% BTBinvsqrt %*% svd..$u
    s <- svd..$d
    tAYA <- t(A) %*% Y %*% A
    yTy <- sum(Y*Y)
    gcv <- numeric(nlam)
    names(gcv) <- format(lambda, digits = 5)
    for (l in 1:nlam) {
        s.tilde <- 1 / (1 + lambda[l]*s)
        s.tilde.sqrt <- sqrt(s.tilde)
        M1 <- outer(s.tilde.sqrt, s.tilde.sqrt) * tAYA
        fran <- sum(M1 * M1)
        M2 <- outer(s.tilde, s.tilde) * tAYA
        joe <- sum(M2 * M2)
        frob2 <- yTy - 2*fran + joe
        gcv[l] <- frob2 / (1 - sum(s.tilde)^2 / n^2)^2
    }
    bestlam <- lambda[which.min(gcv)]
    if (length(lambda)>1 && bestlam %in% range(lambda)) warning("GCV minimized at endpoint of initial grid")
    BBPB <- solve(BTB + bestlam * P, t(B))
    coefmat <- BBPB %*% Y %*% t(BBPB)
    Yhat <- B %*% coefmat %*% t(B)
    list(gcv=gcv, bestlam=bestlam, coefmat=coefmat, Yhat=Yhat)
}


#' Plot an mds.ct object
#' 
#' Plots a curve representing of two continuous-time principal coordinates 
#' produced by \code{\link{mds.ct}}.
#'
#' @param x an object of class "\code{\link{mds.ct}}"
#' @param npoints number of time points (equally spaced along the range of times) 
#' at which to plot the coordinates
#' @param cols color scheme; viridis(50) by default
#' @param title plot title
#' @param size.axes size of axis titles
#' @param samescale logical: Should the coordinates be plotted on the same scale?
#' @param coords which two principal coordinates to plot (default is 1:2)
#' @param GOF.method method to use for computing percent dissimilarity explained
#'  (see argument \code{GOF} of \code{\link{cmdscale}})
#' @param digits number of digits to display for percent dissimilarity explained
#' @param ... other arguments, passed to \code{\link[ggplot2]{theme}}
#'
#' @return None; just produces a plot.
#' 
#' @method plot mds.ct
#' @importFrom fda eval.fd
#' @importFrom ggplot2 ggplot aes theme element_text labs geom_point scale_color_gradientn xlim ylim
#' @importFrom viridisLite viridis
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # see example for mds.ct
#' 
plot.mds.ct <- function(x, npoints=500, cols=viridis(50), title="", size.axes=11,
                        samescale=TRUE, coords=1:2, GOF.method=1, digits=1, ...) {
    if (length(coords) != 2) stop("Must choose two principal coordinates to display")
    pctexp <- 100 * diff(c(0, x$GOF[ , GOF.method]))
    if (is.null(cols)) cols=viridis(50)
    rng <- x$funcs$basis$range
    time <- seq(rng[1], rng[2], , npoints)
    funcvals <- eval.fd(time, x$funcs)[ , coords]
    df <- data.frame(x=funcvals[,1], y=funcvals[,2])
    p <- ggplot(df, aes(x=.data$x, y=.data$y, color=time)) + geom_point() + scale_color_gradientn(colors=cols) +
        labs(x=paste0("Principal coordinate ", coords[1], " (", round(pctexp[coords[1]], digits), "%)"),
             y=paste0("Principal coordinate ", coords[2], " (", round(pctexp[coords[2]], digits), "%)"), 
             title=title) +
        theme(axis.title = element_text(size=size.axes), ...)     
    if (samescale) {
        mn <- min(funcvals); mx <- max(funcvals)
        p <- p + xlim(mn, mx) + ylim(mn, mx)
    }
    p
}


#' Procrustes transformation for continuous-time multidimensional scaling
#' 
#' Matches classical principal coordinates
#' to continuous-time principal coordinates produced by \code{\link{mds.ct}},
#' via Procrustes transformation.
#'
#' @param obj an object of class "\code{\link{mds.ct}}"
#' @param points matrix of classical principal coordinates, 
#' e.g. as produced by \code{\link{cmdscale}}
#' @param coord which coordinates to transform. 
#' 
#' @note
#' The function uses \code{\link[vegan]{procrustes}}, from the \code{vegan}
#' package, to transform the classical 
#' principal coordinates given by \code{points}
#' to the continuous-time principal coordinates defined by \code{obj}, evaluated at the
#' time points given therein. By default, all of the coordinates are used. If \code{obj} and 
#' \code{points} do not include the same number of coordinates, the smaller number
#' is used, and a warning is issued. 
#'
#' @return The output of \code{\link[vegan]{procrustes}}.
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il> 
#'
#' @seealso \code{\link{procrustes.plot.mds}} 
#' @importFrom fda eval.fd
#' @importFrom vegan procrustes
#' @export
#'
#' @examples
#' # see call to plot.procrustes.ct in mds.ct example
#' 
procrustes.mds.ct <- function(obj, points, coord=NULL) {
    argvals <- obj$argvals
    funcs <- obj$funcs
    newpts.all <- eval.fd(argvals, funcs)
    if (is.null(coord)) {
        ncol0 <- ncol(points)
        ncol1 <- ncol(newpts.all)
        if (ncol0 != ncol1) 
            warning("The two MDS representations are of different dimensions")
        coord <- 1:min(ncol0,ncol1)
    }
    procrustes(newpts.all[ , coord], points[ ,coord], scale=FALSE)
}

#' Plot of classical and continuous-time principal coordinates
#' 
#' This function plots two continuous-time principal coordinates along with the 
#' corresponding classical principal coordinates, where the latter is aligned with
#' the former by Procrustes transformation
#'
#' @param obj an object of class "\code{\link{mds.ct}}"
#' @param points matrix of classical principal coordinates, 
#' e.g. as produced by \code{\link{cmdscale}}
#' @param npoints number of time points (equally spaced along the range of times) 
#' at which to plot the coordinates
#' @param cols color scheme; viridis(50) by default
#' @param procoord coordinates to which Procrustes transformation should be applied
#' @param plotcoord which coordinates to plot
#' @param linewidth line width for the principal coordinate curve
#' @param ltitle legend title
#' @param title title of the graph
#' @param samescale logical: Should the x- and y-axes have the same range?
#' @param cex.legend scaling factor for legend key
#' @param GOF.method method to use for computing percent dissimilarity explained
#'  (see argument \code{GOF} of \code{\link{cmdscale}})
#' @param digits number of digits to display for percent dissimilarity explained
#' @param xlim,ylim x- and y-limits. Ignored if \code{samescale==TRUE}.
#' @param size.axis.title size for axis titles 
#' @param size.axis.text size for axis text
#' @param size.ltitle size for color legend title
#' @param size.ltext size for color legend text
#' @param coord_fixed logical: Should aspect ratio be fixed to 1?
#' @param xlab,ylab x- and y-labels. If these are NULL, the principal coordinate numbers 
#' and the percent dissimilarity explained are used as the axis labels.
#' @param colourbar logical: Should a color bar (legend) be included?
#' @param ... arguments passed to \code{\link[ggplot2]{labs}}
#'
#' @return
#' None; a plot is generated
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il> 
#'
#' @seealso \code{\link{procrustes.mds.ct}} 
#' @importFrom fda eval.fd
#' @importFrom ggplot2 geom_path unit labs
#' @importFrom viridisLite viridis
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # see example for mds.ct
#' 
procrustes.plot.mds <- function(obj, points, npoints=500, cols=viridis(50), 
                                procoord=1:2, plotcoord=1:2, linewidth=1.3, ltitle="time", title="",
                                samescale=FALSE, cex.legend=1, GOF.method=1, digits=1,
                                xlim=NULL, ylim=NULL, size.axis.title=11, size.axis.text=11, 
                                size.ltitle=11, size.ltext=11, coord_fixed=TRUE, xlab=NULL, ylab=NULL,
                                colourbar=TRUE, ...) {
    crust <- procrustes.mds.ct(obj=obj, points=points, coord=procoord)
    oldpts <- as.data.frame(matrix(crust$translation[ , plotcoord], nrow(points), 2, byrow=TRUE) +
                                as.matrix(points[ , procoord]) %*% crust$rotation[ , plotcoord])
    newtimes <- seq(min(obj$argvals), max(obj$argvals), , npoints)
    newpts <- as.data.frame(eval.fd(newtimes, obj$funcs)[,plotcoord])
    names(oldpts) <- names(newpts) <- c("x","y")
    oldpts$time <- obj$argvals
    newpts$time <- newtimes
    pctexp <- 100 * diff(c(0, obj$GOF[ , GOF.method]))
    if (is.null(xlab)) 
        xlab <- paste0("Principal coordinate ", plotcoord[1], " (", round(pctexp[plotcoord[1]], digits), "%)")
    if (is.null(ylab)) 
        ylab <- paste0("Principal coordinate ", plotcoord[2], " (", round(pctexp[plotcoord[2]], digits), "%)")
    
    p <- ggplot(oldpts, aes(.data$x,.data$y, color=.data$time)) + geom_point() +
        geom_path(data=newpts, aes(.data$x,.data$y, color=.data$time), linewidth=linewidth) + 
        scale_color_gradientn(colors=cols, name=ltitle, 
                              guide=if (colourbar) "colourbar" else "none") + 
        theme(legend.key.size = unit(cex.legend, "lines"), 
              legend.text  = element_text(size = size.ltext),
              legend.title = element_text(size = size.ltitle),
              axis.title = element_text(size=size.axis.title),
              axis.text = element_text(size=size.axis.text)) + 
        labs(x=xlab, y=ylab, title=title, ...)
    
    if (samescale) {
        if (!is.null(xlim) || !is.null(ylim)) 
            warning("Cannot set 'xlim' or 'ylim' if 'samescale' is TRUE")
        rng <- range(rbind(oldpts,newpts)[,-3])
        p <- p + xlim(rng[1],rng[2]) + ylim(rng[1],rng[2])
    }
    if (!is.null(xlim)) p <- p + xlim(xlim[1], xlim[2])
    if (!is.null(ylim)) p <- p + ylim(ylim[1], ylim[2])
    if (coord_fixed==TRUE) p <- p + coord_fixed(ratio=1)
    p
}

