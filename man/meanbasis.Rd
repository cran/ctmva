% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanbasis.R
\name{meanbasis}
\alias{meanbasis}
\title{Compute means of basis functions}
\usage{
meanbasis(basis, rng = NULL)
}
\arguments{
\item{basis}{a basis object of class \code{"\link[fda]{basisfd}"}}

\item{rng}{time range. By default, the entire interval spanned by the basis. Must be left NULL for Fourier bases.}
}
\value{
Vector of means of the basis functions
}
\description{
Given a basis object as defined in the \pkg{fda} package (see \code{\link[fda]{basisfd}}),
this function simply computes the vector of means of the basis functions. Used internally.
}
\examples{


require(fda)
bbasis6 <- create.bspline.basis(nbasis=6)
meanbasis(bbasis6)
meanbasis(bbasis6, c(.3,.6))
fbasis11 <- create.fourier.basis(nbasis=11)
meanbasis(fbasis11)

}
\author{
Biplab Paul <paul.biplab497@gmail.com> and Philip Tzvi Reiss <reiss@stat.haifa.ac.il>
}
