#' Silhouettes for continuous-time k-means clustering
#' 
#' Computes the silhouette index, at a grid of time points, for a continuous-time 
#' k-means clustering object produced by \code{\link{kmeans.ct}}.
#' 
#' @note An error is issued if the grid of time points contains one or more of the 
#' cluster transition points. This should not ordinarily occur, but if it does, it can 
#' be remedied by modifying \code{ngrid}.
#' 
#' @param kmobj continuous-time k-means clustering from \code{\link{kmeans.ct}}
#' @param ngrid number of equally spaced grid points at which to compute the silhouette index
#' @return Object of class "\code{silhouette.ct}", a list consisting of
#' \item{grid }{grid of \code{ngrid} points spanning the time range} 
#' \item{value }{silhouette index at each point along the grid} 
#' \item{transitions }{transition points between segments}
#' \item{cluster }{cluster memberships in the segments defined by the transitions}
#' \item{mean }{mean silhouette index}
#' 
#' @author Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
#' @seealso \code{\link{kmeans.ct}}, which includes an example; \code{\link{plot.silhouette.ct}}
#' @export silhouette.ct
silhouette.ct <- function(kmobj, ngrid=5000) {
    k <- max(kmobj$cluster)
    rng <- kmobj$fdobj$basis$rangeval
    starts <- c(rng[1], kmobj$transitions)
    ends <- c(kmobj$transitions, rng[2])
    grid <- seq(rng[1], rng[2], , ngrid)
    grid2 <- sort(c(grid, kmobj$transitions))
    C <- kmobj$fdobj$coef
    Phi <- fda::eval.basis(grid2, kmobj$fdobj$basis)
    D <- as.matrix(stats::dist(Phi %*% C))[which(grid2 %in% grid), ]
    if (nrow(D) != ngrid) stop("grid includes transition points; please modify ngrid")
    
    ai <- bi <- cluster_long <- rep(0, ngrid)
    md2c <- matrix(NA, ngrid, k)  # mean distance to each cluster

    for (m in 1:k) {
        
        which_sections <- which(kmobj$cluster==m)

        # For each t in grid2, integrate distance-from-t over cluster-m sections
        int_d <- 0
        for (section in which_sections) {
            cluster_long[grid >= starts[section] & grid <= ends[section]] <- m
            
            # Integrate distance-from-t within section
            which_pts <- grid2 >= starts[section] & grid2 <= ends[section]
            pts <- grid2[which_pts]
            delta <- diff(pts)
            wts <- (c(0,delta) + c(delta,0))/2  # trapezoidal rule weights
            int_section <- apply(D, 1, function(v) sum(v[which_pts] * wts))
            int_d <- int_d + int_section
        }
        
        md2c[ , m] <- int_d / kmobj$size[m]
    }
    
    abc <- md2c
    for (m in 1:k) {
        which_m <- (cluster_long==m)
        ai[which_m] <- md2c[which_m, m]
        abc[which_m, m] <- max(md2c) + 1
    }
    bi <- apply(abc,1,min)
    lst <- list(grid=grid, value=(bi-ai)/pmax(ai,bi), 
                transitions=kmobj$transitions, cluster=kmobj$cluster)
    lst$mean <- mean(lst$value)
    class(lst) <- "silhouette.ct"
    lst
}
