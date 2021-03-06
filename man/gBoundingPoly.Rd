% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mSpatial.R
\name{gBoundingPoly}
\alias{gBoundingPoly}
\title{bounding box polygon function}
\usage{
gBoundingPoly(container, tolerance = 1.5)
}
\arguments{
\item{container}{The single polygon to whose extent the reference polygons will be included}

\item{tolerance}{The multiple by which the container extent will be expanded (default = 1.5)}
}
\description{
This function returns a bounding box polygon whose extent is some multiple of the extent of a container polygon. It's most useful as an input to other reducing operations.
}
\details{
\strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
}
\examples{
gBoundingPoly(drivetime, tolerance = 1.25)
}
\keyword{spatial}
