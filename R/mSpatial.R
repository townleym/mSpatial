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
gBoundingPoly = function(container, tolerance = 1.5) {

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
			
		bp = as(raster::extent(newext), "SpatialPolygons")	
		proj4string(bp) = proj4string(container)
		bp
	}
}

#' bounding box clip function
#' 
#' This function clips a set of reference polygons to some multiple of the extent of a container polygon. It is most useful for reducing a large set of polygons (e.g. US Census Block groups) to a smaller subset. \strong{New!} it now returns full intersecting polygons rather than the clipped slivers. And if reference is a SpatialPolygonsDataFrame it will return a SpatialPolygonsDataFrame (yay!)
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
        proj4string(temp_bp) = CRS(proj4string(container))		
		reference[raster::intersect(reference, temp_bp),]
		# rgeos::gIntersection(reference, temp_bp, byid = T)
	}
}

#' Return polygons enclosed by a container --based on a percent area enclosed (deprecated)
#' 
#' We want to return the reference polygons whose areal overlap with the container falls within some threshold
#'
#' This function uses name matching / extraction. It's a kludge and often fails. use \code{gPolyByIntersect} instead. Preserved for existing code that uses it.
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

#' Return polygons enclosed within a container (deprecated)
#' 
#' Pretty much a copy of \code{gEnclosedCentroid} but more robust
#'
#' Deprecated in favor \code{gPolyByIntersect(container, reference, centroid = T)} instead. Preserved for existing code that uses it.
#'
#' Returns reference polygons whose centroids are enclosed within a containing polygon. Allows for (requires) the pre-specification of reference polygon centroids which is useful when there are a very large number of reference polygons. Uses a bounding box filter before any spatial operations which makes it much faster than a straight intersection.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param reference A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} object with the values to be summarized within a containing polygon
#' @param centroids A set of centroids for the reference polygons obtained by \code{gCentroids}
#' @param container A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} object of a single, container polygon
#' @keywords spatial
#' @export
#' @examples 
#' gPolyByCentroid(blockgroups, sp::gCentroid(blockgroups, byid = T), cbsa)
gPolyByCentroid = function(reference, centroids, container) {

	# 1 Use centroids of incoming polys to cut down to a bounding box
	ta_blocks_bbox_centroids_idx = gWhichPoints(centroids, bbox(container))
	# 2 Cut incoming polys down to those in the bounding box
	ta_blocks_bbox = reference[ta_blocks_bbox_centroids_idx,]
	# 3 return polys (from the bbox above) whose centroids lie in the container
	ta_blocks_bbox[which(gContains(container, gCentroid(ta_blocks_bbox, byid = T), byid = T)),]

}


#' Error Trap
#' 
#' Copied from Hadley. Useful for providing a way for functions to exit gracefully from an error
#'
#' Returns \code{TRUE} if \code{try} returns an error
#'
#' @param any value returned from a function
#' @keywords none
#' @export
#' @examples 
#' result = try(somefunction())
#' if(is.error(result)) {
#'	print("Bailing out!")
#' } else {
#'	print("result is not an error!")
#'	}
is.error = function(x) {
	inherits(x, "try-error")
}

#' Return polygons enclosed within a container (Deprecated)
#' 
#' Kinda like \code{gEnclosedArea}
#'
#' Returns all the reference polygons from a SpatialPolygonsDataFrame that intersect with a container. Adds an attribute to the @data slot with the percent of overlapping area. 
#'
#' This function uses name matching / extraction. It's a kludge and often fails. use \code{gPolyByIntersect} instead. Preserved for existing code that uses it.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|Points}
#;
#' @param reference A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} object with the values to be summarized within a containing polygon
#' @param centroids A set of centroids for the reference polygons obtained by \code{gCentroids}
#' @param container A \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} object of a single, container polygon
#' @keywords spatial
#' @export
#' @examples 
#' gPolyByCentroid(blockgroups, sp::gCentroid(blockgroups, byid = T), cbsa)

# This time with error handling!
# maybe give option to return slivers or full overlap
mGetIntersecting = function(container, reference) {

	clipped = try(rgeos::gIntersection(reference, container, byid = T))

	if (is.error(clipped) | is.null(clipped)) { # catch errors in the gIntersection() function
		return(NA)
	} else {
		clipped_names = sapply(strsplit(names(clipped), " "), "[[", 1)
	
		# calc area of intersecting parts
		clipped_area = clipped %>% gArea(byid = T)
		names(clipped_area) = sapply(strsplit(names(clipped), " "), "[[", 1)
		clipped_area = data.frame(clipped_area = clipped_area, container_overlap = 	clipped_area / gArea(container))
		clipped_reference = reference[which(rownames(reference@data) %in% clipped_names),]
		
		retframe = sp::merge(clipped_reference, clipped_area, by = "row.names")
		return( retframe[order(retframe$container_overlap, decreasing = T),] )
	} # end exception trap
} # end mGetIntersecting

#' Extract X, Y coordinates from a WKT points object
#'
#' Quick little thingy to run over a list of WKT points and convert to a data frame of lat/lon. Mostly done because some people think it's really clever to.... yeah
#;
#' @param ptwkt A bunch of WKT formatted single points
#' @keywords spatial
#' @export
#' @examples 
#' mWKTPointToXY(list_of_WKT_points)

# this seems unnecessarily tedious
mWKTPointToXY = function(ptwkt) {
    thing = lapply(ptwkt, function(x) coordinates(readWKT(x)))
    exes = sapply(thing, "[[", 1)
    whys = sapply(thing, "[[", 2)
    data.frame(lon = exes, lat = whys)
}

#' Return a spatial polygon from a point
#' 
#' Returns an object of type \code{SpatialPolygons} with a single feature. 
#' 
#' Default behavior is a square polygon based on radial distance (specified by \code{size}) in the unit of the projection spatial reference system. Default linear measure is in meters. Because it is defined radially a square buffer will have sides equal to 2x the \code{size} argument. 
#' 
#' If \code{square = F} the \code{size} argument will be scaled as a multiple of \code{long.ratio}. If \code{landscape = T} (default) is chosen, then the X-axis (horizontal) side will be longer. Note that \code{square = F} implies a 1:1 correspondence between \code{size} and the shortest side of the rectangle. In other words, if the \code{size} argument is unchanged, the short side of a box from \code{square = F} will be one half the length of the side of a box where \code{square = T}. I could fix that pretty easily. But this works fine for me...
#' 
#' The function will return a polygon in the same coordinate system as the point given as input (specified as \code{gcs.srid}). It makes no attempt to transform the spatial reference system of the input point. As of now, it's only been tested with points in WGS84 (\code{epsg:4326}) lat/long (x = long, y = lat).
#' 
#' \strong{Note:} The default (flawed) "Web Mercator" (\code{EPSG:3857}) was chosen as the planar spatial reference system so as not to bias in favor of any world region. Since it is used here for area calculations and given the Mercator's well-known exaggeration of areas as latitudes get larger, an equal-area projection is recomended (e.g. in the U.S., the US National Atlas \code{"+init=epsg:2163"}).
#'
#' @param x Real value corresponding to longitude of a point (expected as WGS84 or \code{epsg:4326})
#' @param y Real value corresponding to latitude of a point (expected as WGS84 or \code{epsg:4326})
#' @param wktstring a WKT formatted point \strong{required} if \code{x, y} not given
#' @param size radial size of box given in units for the given projection (see Details)
#' @param square should a square box (of area 2*size^2) be returned; if \code{F} then a rectangle whose long side is a multiple of size configured by...
#' @param long.ratio The factor for stretching the box along its long axis, configured by....
#' @param landscape whether to make the long dimension the horizontal (x) (default) or vertical
#' @param proj.srid a \code{proj4string} for a planar projected coordinate system. Default \code{epsg:3857} or "Web Mercator"
# #' @param gcs.srid a \code{proj4string} for the input geographic (lat/lon) coordinate system and for the returned object. Default \code{epsg:4326} or "WGS84"
#' @keywords spatial
#' @export
#' @examples 
#' box_1k = mBox(wktstring = "POINT(90 30)", size = 500)
#' box_1k = mBox(-90, 30, size = 500)
#' rectangle_1k = mBox(-90, 30, size = 1e3, square = F, long.ratio = 16/9)
mBox = function(x, y, wktstring, size = 500, square = T, long.ratio = 1.618034, landscape = T, proj.srid = "+init=epsg:3857", gcs.srid = "+init=epsg:4326") {
		
	if(!exists("wktstring")) {
		wktstring = paste0("POINT(", x, " ", y,")")	
	}
	# convert the x,y coordinates to a spatial point with a useful 
	# planar projection
	pint = readWKT(wktstring, p4s = CRS(gcs.srid))
	pint_proj = spTransform(pint, CRSobj = CRS(proj.srid))

	# Calculate the extent of the box either as a square...
	if(square) {
		ext = gBuffer(pint_proj, width = size) %>% bbox
		
	} else {	# ...or as a rectangle with one side longer according to long.ratio
		long.side = size * long.ratio
	
		x.coord = coordinates(pint_proj)[,"x"]
		y.coord = coordinates(pint_proj)[,"y"]

		if(landscape) {
			exes = x.coord + ((c(-1,1) * (long.side / 2)))
			whys = y.coord + ((c(-1,1) * (size / 2))) 
		} else {
			exes = x.coord + ((c(-1,1) * (size / 2)))
			whys = y.coord + ((c(-1,1) * (long.side / 2))) 
		}
		ext = as.matrix(rbind(x = exes, y = whys))
		colnames(ext) = c("min", "max")

	} # end rectangle

	prj_box = as(raster::extent(ext), "SpatialPolygons")
	proj4string(prj_box) = proj.srid
	spTransform(prj_box, CRSobj = CRS(gcs.srid))

} # end mBox function


#' Convert R color strings to hex 
#' 
#' Quickly create a hex string for a color + alpha channel
#' 
#' @param basecol An R base color name
#' @param achannel (ideally) a hex value string for the alpha channel. Decimal values will be interpreted as hex
#' 
#' @export
#' @examples
#' col2hex("lightsteelblue", "70")
#' sapply(brewer.pal(n = 5, name = "BuPu"), function(x) {col2hex(x, achannel = '80')})
col2hex = function(basecol, achannel = "") {
	col2rgb(basecol) %>% as.hexmode %>% paste0(collapse = "") %>% paste0("#", ., achannel, collapse = "")
}

#' Return polygons intersecting with a container
#' 
#' This is the preferred function. It has been optimized for speed and flexibility. Allows choice of centroid or areal overlap (intersection). Amount of overlap (for non-centroid intersection methods) specified by \code{threshold} defaults to zero which returns an intersection for two polygons that share any one point (unfortunately) including a border.
#'
#' \strong{Note}: all inputs should be of class \code{sp::SpatialPolygons|SpatialPolygonsDataFrame}.
#' 
#' Will warn if CRS mismatch between reference / container. 
#' 
#' @param container The single polygon that will enclose the set of reference polygons
#' @param reference The set of polygons that will be selected down to the containing polygon based on areal overlap
#' @param threshold The percent of areal overlap required for inclusion (ignored if \code{centroid = T})
#' @param centroid Logical return reference polygons whose centroids are contained within the container
#' @keywords spatial
#' @export
#' @examples 
#' bgoi.clip = gBoundingPolyClip(box_parkway, bg_no_data)
#' bgoi.cent = gPolyByIntersect(box_parkway, bgoi.clip, centroid = T)
#' bgoi.area.20 = gPolyByIntersect(box_parkway, bgoi.clip, threshold = 0.2) 
#' bgoi.intersect = gPolyByIntersect(box_parkway, bgoi.clip) 
gPolyByIntersect = function(container, reference, threshold = 0, centroid = F) {
	if(centroid) {
		
		centroids = gCentroid(reference, byid = T)
		# 1 Use centroids of incoming polys to cut down to a bounding box
		ta_blocks_bbox_centroids_idx = gWhichPoints(centroids, bbox(container))
		# 2 Cut incoming polys down to those in the bounding box
		ta_blocks_bbox = reference[ta_blocks_bbox_centroids_idx,]
		# 3 return polys (from the bbox above) whose centroids lie in the container
		ta_blocks_bbox[which(gContains(container, gCentroid(ta_blocks_bbox, byid = T), byid = T)),]
		# end byCentroid

	} else {
		
		slivers = raster::intersect(reference, container)
		# use the google web mercator projection for area calc
		sliver.area = gArea(spTransform(slivers, CRS("+init=epsg:3857")), byid = T) 
		reference.area = gArea(spTransform(reference[slivers,], CRS("+init=epsg:3857")), byid = T)
		overlap = sliver.area / reference.area
		reference[slivers,][which(overlap > threshold),]

	} # end byOverlap
} # end gPolyByIntersect

#' Plot OSM basemap (deprecated)
#' 
#' Convenience wrapper around the openmap() function to return the basemap around an input geometry. This has been deprecated in favor of \code{topleft()} and \code{bottomright()}.
#'
#' \strong{Note}: for now it really only works around a polygon
#' 
#' @param geom Input geometry (polygon)
#' @param tolerance multiplier for the size of the bounding box around the input geometry
#' @param tileserver see the openmap() documentation for options
#'
#' @keywords spatial
#' @export
#' @examples 
#' basemap = osmGet(cambridge_town, tolerance = 1.5, tileserver = "mapbox")
#' plot(basemap)
#' # skobbler maps look better with a little purple haze
#' plot(spTransform(gBoundingPoly(cambridge_town, tolerance = 1.5), osm()), add = T, col = col2hex("lavender", "40"), border = NA)
#' plot(spTransform(cambridge_town, osm()), add = T, col = col2hex("skyblue", "80"))
osmGet = function(geom, tolerance = 2, tileserver = "skobbler") {
		
	if(require(OpenStreetMap)) {
		ext_map_box = gBoundingPoly(geom, tolerance)
		# proj4string(ext_map_box) = CRS("+init=epsg:4326")
		proj4string(ext_map_box) = CRS(proj4string(geom))
		basemap_ext = bbox(ext_map_box)

		basemap_topleft = basemap_ext %>% diag %>% rev 
		basemap_bottomright = c(basemap_ext[2,1], basemap_ext[1,2])

		openmap(basemap_topleft, basemap_bottomright, type = tileserver)	
	}	
}


#' Retrieve topleft of \code{sp} object
#' 
#' For use with the (annoying) \code{openmap()} functions which retrieve basemap tiles for a box specified by the topleft and bottom right vertices.
#'
#' @param geom Input geometry (polygon)
#'
#' @keywords spatial
#' @export
#' @examples 
#' shade = gBoundingPoly(cambridge_town, tol = 1.1)
#' basemap = openmap(topleft(shade), bottomright(shade), type = "skobbler")
#' 
#' plot(basemap)
#' # skobbler maps look better with a little purple haze
#' plot(spTransform(shade, osm()), add = T, col = col2hex("lavender", "40"), border = NA)
#' plot(spTransform(cambridge_town, osm()), add = T, col = col2hex("skyblue", "80"))
topleft = function(geom) {
	box = bbox(geom)
	c(max(box[2,]), min(box[1,]))	
}

#' Retrieve bottomright of \code{sp} object
#' 
#' For use with the (annoying) \code{openmap()} functions which retrieve basemap tiles for a box specified by the topleft and bottom right vertices.
#'
#' @param geom Input geometry (polygon)
#'
#' @keywords spatial
#' @export
#' @examples 
#' shade = gBoundingPoly(cambridge_town, tol = 1.1)
#' basemap = openmap(topleft(shade), bottomright(shade), type = "skobbler")
#' 
#' plot(basemap)
#' # skobbler maps look better with a little purple haze
#' plot(spTransform(shade, osm()), add = T, col = col2hex("lavender", "40"), border = NA)
#' plot(spTransform(cambridge_town, osm()), add = T, col = col2hex("skyblue", "80"))
bottomright = function(geom) {
	box = bbox(geom)
	c(min(box[2,]), max(box[1,]))
}

#' summarize attributes over a containing geometry
#' 
#' Convenience wrapper around gPolyByIntersect that applies a function over specified attributes in the reference set. It allows for buffering the input container.
#'
#' \strong{Note}: careful, it's pretty fragile
#' 
#' @param container The thing you want to summarize within
#' @param reference The thing whose attributes you want to summarize
#' @param buffer How much you want to blow up the summarizing container (in meters)
#' @param colnamevec A vector of column names corresponding to those in the container whose values will be summarized
#' @param centroid Whether to use containing centroids for choosing reference geometries within the container
#' @param threshold the overlap percentage to use when \code{centroid = F}
#'
#' @keywords spatial
#' @export
#' @examples 
#' mSummarizer(cambridge_town, boston_blocks, buffer = 0, c("total_pop", "total_hh"), centroid = T, FUN = "sum")
mSummarizer = function(container, reference, buffer = 1, colnamevec, centroid = F, threshold = 0.99, FUN) {
	
	# constant
	gmap = "+init=epsg:3857"
	# set incoming container projection to back transform after projection	
	container_proj = proj4string(container)

	# for dev
	# container = beacon_mi
	# reference = bg_hh_car
	# buffer = 1
	# colnamevec = names(hh_car)[c(3,5,6,8,9)]
	# centroid = T
	# threshold = 0.99
	# FUN = "sum"
	
	# Get the right function to apply()
	FUN = match.fun(FUN)
	
	# buffer the incoming geography (defaults to no buffer)
	buffered_container = gBuffer(spTransform(container, CRSobj = gmap), width = buffer) %>% spTransform(CRSobj = CRS(container_proj))
	
	# Get the reference polygons contained within the container
	sub_ref = gPolyByIntersect(buffered_container, reference, centroid = centroid, threshold = threshold)

	# Sum up attributes given by colnamevec	
	apply(sub_ref@data[,colnamevec], 2, FUN)
}

#' Height / width ratio
#' 
#' Quick function for getting the 1/aspect ratio of a geometry
#' 
#' @param geom an \code{sp} geometry object
#' 
#' Since map objects have varying aspect ratios, saving them to files results in lots of whitespace. This function enables you to specify the output dimensions in the same ratio as the \code{sp} object.
#' 
#' @export
#' @examples
#' x = 800 # width in pixels
#' y = yoverx(basemap)
#' png(filename, width = x, height = y)
yoverx = function(geom, osm = T) {
	
	if(osm) {
		ta.bbox = bbox(spTransform(geom, osm()))
		xy = ta.bbox[,"max"] - ta.bbox[,"min"]
		xy[2] / xy[1]		
		
	} else {
		laea = "+init=epsg:2163" # laea US natl atlas
	
		ta.bbox = bbox(spTransform(geom, CRSobj = CRS(laea)))
		xy = ta.bbox[,"max"] - ta.bbox[,"min"]
		xy[2] / xy[1]		
	}
}
