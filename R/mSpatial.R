#' bounding box polygon function
#' 
#' This function returns a bounding box polygon whose extent is some multiple of the extent of a container polygon. It's most useful as an input to other reducing operations.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param container The single polygon to whose extent the reference polygons will be included
#' @param tolerance The multiple by which the container extent will be expanded (default = 1.5)
#' @keywords spatial
#' @export
#' @examples 
#' gBoundingPoly(drivetime, tolerance = 1.25)

gBoundingPoly = function(container, reference, tolerance = 1.5) {

	if(require(raster)) {
		expansion.factor = (tolerance - 1) / 2
	
		# get the extent of the container polygon	
		v.ext = as.vector(t(sp::bbox(container)))
		names(v.ext) = c("xmin", "xmax", "ymin", "ymax")
	
		xrange = v.ext["xmax"] - v.ext["xmin"]
		yrange = v.ext["ymax"] - v.ext["ymin"]
	
		# Calculate the new extent	
		newext = bbox(container)
		newext[1,] = v.ext[1:2] + c(-1,1) * (expansion.factor * xrange)
		newext[2,] = v.ext[3:4] + c(-1,1) * (expansion.factor * yrange)
			
		as(raster::extent(newext), "SpatialPolygons")	
	}
}

#' bounding box clip function
#' 
#' This function clips a set of reference polygons to some multiple of the extent of a container polygon. It is most useful for reducing a large set of polygons (e.g. US Census Block groups) to a smaller subset.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param container The single polygon to whose extent the reference polygons will be included
#' @param reference The set of polygons that will be selected down to the reference polygon
#' @param tolerance The multiple by which the extent of the reference polygon will be expanded (default = 1.5)
#' @keywords spatial
#' @export
#' @examples 
#' gBoundingPolyClip(drivetime, block_groups, tolerance = 1.25)

gBoundingPolyClip = function(container, reference, tolerance = 1.5) {

	if(require(raster)) {
		# container = roggies_drivetime
		# reference = ma_bg
		# tolerance = 2

		# clip reference layer to some extent related to the container polygon
		expansion.factor = (tolerance - 1) * (1/2)
	
		# get the extent of the container polygon	
		ext = sp::bbox(raster::extent(container))
		xrange = ext[1,2] - ext[1,1]
		yrange = ext[2,2] - ext[2,1]
	
		# Calculate the new extent	
		newext = ext
		newext[1,] = ext[1,] + c(-1,1) * (expansion.factor * xrange)
		newext[2,] = ext[2,] + c(-1,1) * (expansion.factor * yrange)
			
		temp_bp = as(raster::extent(newext), "SpatialPolygons")	
		rgeos::gIntersection(reference, temp_bp, byid = T)
	}
}

#' Return polygons enclosed by a container --based on a percent area enclosed
#' 
#' We want to return the reference polygons whose areal overlap with the container falls within some threshold
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}.
#' 
#' Areal calculations really only work if the input polygons are rendered in the same planar projection.
#'
#' @param container The single polygon that will enclose the set of reference polygons
#' @param reference The set of polygons that will be selected down to the containing polygon based on areal overlap
#' @param threshold The percent of areal overlap required for inclusion
#' @keywords spatial
#' @export
#' @examples 
#' require(sp)
#' # Use the US National Atlas Lambert Equal Area projection
#' drivetime.proj = spTransform(drivetime, CRSobj = CRS(EPSG[EPSG$code == 2163 & !is.na(EPSG$code), "prj4"]))
#' blockgroups.proj = spTransform(block_groups, CRSobj = CRS(EPSG[EPSG$code == 2163 & !is.na(EPSG$code), "prj4"]))
#' # Return only the blockgroups whose overlap exceeds 20%
#' gEnclosedArea(drivetime.proj, block_groups.proj, threshold = 0.2)
#' 
#' # Often it will be faster to work with a subset of reference polygons 
#' reference.bg.proj = gBoundingPoly(drivetime.proj, blockgroups.proj, tolerance = 1.25)
#' gEnclosedArea(drivetime.proj, reference.bg.proj, threshold = 0.2)
gEnclosedArea = function(container, reference, threshold = 0.99) {

	# reference is a set of whole polygons for which you have some attributes
	# clipped are the reference polygons 'clipped' to the container (not used in this function)
	# 
	# We want to return the reference polygons whose overlap with the container
	# falls within some threshold
	# 
	# This is all done by matching names of sets of polygons then extracting by index
	# Lots of index/name swapping
	#
	# This function evaluates whether a reference polygon should be in our result set
	# by a threshold of area covered by the larger container
	# 
	# 1. clipped = clip (gIntersection) reference to container
	# 2. Get the names/indices of the reference polygons within the container (clipped)
	# 3. Calculate areas of the contained areas of the clipped polygons
	# 4. Calculate areas of the complete reference polygons (that have a match in the clipped)
	# 5. Find names/indices whose (clipped area / reference area) is > threshold
	# 6. Return the whole polygons from the reference set that meet the threshold
	
	# 1. clip the reference polygons to the container
	clipped = rgeos::gIntersection(reference, container, byid = T)
		
	# 2. match by name
	clipped_bg_names = sapply(strsplit(names(clipped), " "), "[[", 1)
	reference_bg_names = sapply(strsplit(names(reference), " "), "[[", 1)
	
	# matching ids of every reference polygon that has any intersect with the polygons
	# clipped to the container
	all_match_ids = which(reference_bg_names %in% clipped_bg_names)

	# 3 & 4. Get ratio of areas
	areas = round(rgeos::gArea(clipped, byid = T) / rgeos::gArea(reference[all_match_ids,], byid = T), 3)
	
	# 5. Find names then indices of reference polygons that meet threshold
	threshold_match_names = sapply(strsplit(names(areas[areas > threshold]), " "), "[[", 1)		
	threshold_match_ids = which(reference_bg_names %in% threshold_match_names)

	# 6. Return the subset of reference polygons whose intersecting area meets the threshold
	reference[threshold_match_ids,]
}


#' Return polygons whose centroids are enclosed 
#' 
#' We want to return the reference polygons whose centroids are contained within a larger, containing polygon
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param container The single polygon that will enclose the set of reference polygons
#' @param reference The set of polygons that will be selected down to the containing polygon
#' @keywords spatial
#' @export
#' @examples 
#' gEnclosedCentroid(drivetime, block_groups)
gEnclosedCentroid = function(container, reference) {
	something = rgeos::gContains(container, rgeos::gCentroid(reference, byid = T), byid = T)
	# use subsetting to return polygons whose centroids are contained within the above polygon
	reference[something,]
}

#' Select points within bounding box
#' 
#' Uses bracket extraction to very quickly select points within a bounding box
#'
#' Intersect methods (e.g. \code{rgeos::gContains} or \code{rgeos::gIntersection}) are very slow for extracting points via a bounding box. This function simply subsets a \code{SpatialPoints} or \code{SpatialPointsDataFrame} object using standard bracket extraction. As a result, it is much faster
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param allpoints A \code{SpatialPoints} or \code{SpatialPointsDataFrame} object
#' @param b_box A set of bounding coordinates in matrix form like those returned from \code{sp::bbox()}
#' @keywords spatial
#' @export
#' @examples 
#' gBoundingPoints(points.spdf, bbox(container_polygon))
gBoundingPoints = function(allpoints, b_box) {

	if(class(b_box) == "matrix") {
		v.ext = as.vector(t(b_box))
		names(v.ext) = c("xmin", "xmax", "ymin", "ymax")

		allpoints[allpoints@coords[,1] > v.ext["xmin"] & # longitude
		allpoints@coords[,1] < v.ext["xmax"] & # longitude
		allpoints@coords[,2] > v.ext["ymin"] & # latitude
		allpoints@coords[,2] < v.ext["ymax"],] # latitude
	}
	else {
		stop("What kind of bounding box are you throwing at me, smalls?")
	}
}

#' Return indices of points within a box
#' 
#' Like \code{gBoundingPoints} but only returns the indices of points contained within a bounding box
#'
#' Rather than returning the spatial object, this give the indices of points within a bounding box. It's more useful for re-merging attributes (Since gBoundingPoints strips a SpatialPointsDataFrame to a SpatialPoints object). Or, better still, subsetting.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param allpoints A \code{SpatialPoints} or \code{SpatialPointsDataFrame} object
#' @param b_box A set of bounding coordinates in matrix form like those returned from \code{sp::bbox()}
#' @keywords spatial
#' @export
#' @examples 
#' gWhichPoints(points.spdf, bbox(container_polygon))
#' all_the_points_spdf[gWhichPoints(all_the_points.spdf, bbox(container_polygon)), ]
gWhichPoints = function(allpoints, b_box) {

	if(class(b_box) == "matrix") {
		
		# allpoints = ma_lm_point
		# b_box = bbox(boston)

		v.ext = as.vector(t(b_box))
		names(v.ext) = c("xmin", "xmax", "ymin", "ymax")
		which(
		allpoints@coords[,1] > v.ext["xmin"] & # longitude
		allpoints@coords[,1] < v.ext["xmax"] & # longitude
		allpoints@coords[,2] > v.ext["ymin"] & # latitude
		allpoints@coords[,2] < v.ext["ymax"] # latitude
		)
	}
	else {
		stop("What kind of bounding box are you throwing at me, smalls?")
	}
}
