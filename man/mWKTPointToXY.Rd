% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mSpatial.R
\name{mWKTPointToXY}
\alias{mWKTPointToXY}
\title{Extract X, Y coordinates from a WKT points object}
\usage{
mWKTPointToXY(ptwkt)
}
\arguments{
\item{ptwkt}{A bunch of WKT formatted single points}
}
\description{
Quick little thingy to run over a list of WKT points and convert to a data frame of lat/lon. Mostly done because some people think it's really clever to.... yeah
}
\examples{
mWKTPointToXY(list_of_WKT_points)
}
\keyword{spatial}
