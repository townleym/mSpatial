% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mSpatial.R
\name{gPolyByIntersect}
\alias{gPolyByIntersect}
\title{Return polygons intersecting with a container}
\usage{
gPolyByIntersect(container, reference, threshold = 0, centroid = F)
}
\arguments{
\item{container}{The single polygon that will enclose the set of reference polygons}

\item{reference}{The set of polygons that will be selected down to the containing polygon based on areal overlap}

\item{threshold}{The percent of areal overlap required for inclusion (ignored if \code{centroid = T})}

\item{centroid}{Logical return reference polygons whose centroids are contained within the container}
}
\description{
This is the preferred function. It has been optimized for speed and flexibility. Allows choice of centroid or areal overlap (intersection). Amount of overlap (for non-centroid intersection methods) specified by \code{threshold} defaults to zero which returns an intersection for two polygons that share any one point (unfortunately) including a border.
}
\details{
\strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|SpatialPolygonsDataFrame}.

Will warn if CRS mismatch between reference / container.
}
\examples{
bgoi.clip = gBoundingPolyClip(box_parkway, bg_no_data)
bgoi.cent = gPolyByIntersect(box_parkway, bgoi.clip, centroid = T)
bgoi.area.20 = gPolyByIntersect(box_parkway, bgoi.clip, threshold = 0.2) 
bgoi.intersect = gPolyByIntersect(box_parkway, bgoi.clip) 
}
\keyword{spatial}
