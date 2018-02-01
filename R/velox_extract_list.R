#' @title Extract Values Given Polygons
#'
#' @name VeloxRaster_extract_list
#'
#' @description
#' Extracts the values of all cells intersecting with a spatial object (line or polygon)
#' \code{sp} and optionally applies R function \code{fun}.
#'
#' @details
#' If passed, \code{fun} must be an R function accepting a numeric vector as its first (and only mandatory) argument, and returning a anything that can in a list.
#' If \code{fun} is \code{NULL}, \code{extract} returns a list of matrices, each matrix containing the raster values intersecting with the respective polygon (but see argument \code{df}).
#' If sp contains polygons, then cell-polygon intersections are calculated based on cell centroids (but see argument \code{small}).
#' If sp contains lines, then regular cell-line intersections are calculated.
#'
#' @param sp A sf* POLYGON or MULTIPOLYGON object, a sf* LINE or MULTILINE object, a SpatialPolygons* object, or a SpatialLines* object.
#' @param fun An R function. See Details.
#' @param df Boolean. If TRUE, the return value will be a data frame (or list of data frames, see Details), otherwise a matrix (or list of matrices, see Details).
#' If TRUE, a column \code{ID_sp} will be added to each data frame containing the ID of the sp object.
#' @param small Boolean. If TRUE and sp contains polygons, then raster values for small (or oddly shaped) polygons that do not intersect with any cell centroid
#' are established by intersecting the small polygon with the entire (boxed) cells.
#' @param legacy Boolean. Whether to use legacy C++ code (pre velox 0.1.0-9007).
#'
#' @return
#' If \code{fun} is passed: A numeric matrix or data frame (see argument \code{df}) with one row per element in \code{sp}, one column per band in the VeloxRaster.
#'
#' Otherwise: A list of numeric matrices or data frames (see argument \code{df}), with one list element per element in \code{sp}.
#' Each matrix/data frame consists of one column per band in the VeloxRaster, one row per raster cell intersecting with the geometry.
#'
#' @examples
#' ## Make VeloxRaster with two bands
#' set.seed(0)
#' mat1 <- matrix(rnorm(100), 10, 10)
#' mat2 <- matrix(rnorm(100), 10, 10)
#' vx <- velox(list(mat1, mat2), extent=c(0,1,0,1), res=c(0.1,0.1),
#'             crs="+proj=longlat +datum=WGS84 +no_defs")
#' ## Make SpatialPolygons
#' library(sp)
#' library(rgeos)
#' coord <- cbind(0.5, 0.5)
#' spoint <- SpatialPoints(coords=coord)
#' spols <- gBuffer(spgeom=spoint, width=0.5)
#' ## Extract
#' vx$extract(sp=spols, fun=mean)
#'
#' @import rgeos
#' @import sp
NULL
VeloxRaster$methods(extract_list = function(sp, fun = NULL, small = FALSE) {
  "See \\code{\\link{VeloxRaster_extract_list}}."


    # Ensure we have sfc object
    if (inherits(sp, "sf")) {
      sp <- st_geometry(sp)
    }

    # Convert to sfc / ensure argument class is correct
    isLine <- FALSE
    if (inherits(sp, "SpatialPolygons") || inherits(sp, "SpatialPolygonsDataFrame")) {
      geomc <- st_geometry(st_as_sf(sp))
      sp_IDs <- sapply(sp@polygons, function(x) slot(x, "ID"))
    } else if (inherits(sp, "SpatialLines") || inherits(sp, "SpatialLinesDataFrame")) {
      geomc <- st_geometry(st_as_sf(sp))
      sp_IDs <- sapply(sp@lines, function(x) slot(x, "ID"))
      isLine <- TRUE
    } else if (inherits(sp, "sfc_MULTIPOLYGON") || inherits(sp, "sfc_POLYGON")) {
      geomc <- sp
      sp_IDs <- 1:length(geomc)  # No polygon IDs in sf...?
    } else if (inherits(sp, "sfc_MULTILINESTRING") || inherits(sp, "sfc_LINESTRING")) {
      geomc <- sp
      sp_IDs <- 1:length(geomc)  # No polygon IDs in sf...?
      isLine <- TRUE
    } else {
      stop("Argument sp is of wrong class.")
    }

    # Prepare out
    if (!is.null(fun)) {
		# see: https://stackoverflow.com/a/30007964/140384
      out <- matrix(rep(list(),NROW(geomc)*nbands), NROW(geomc), nbands)
      rownames(out) <- sp_IDs
    } else {
      out <- vector("list", NROW(geomc))
      names(out) <- sp_IDs
    }

    # Ensure we have an overlap
    overlaps <- .self$overlapsExtent(geomc)
    if (!overlaps) {
      return(out)
    }

    # Create boost grid, boost geometries, intersect
    if (isLine) {
      boostGrid <- boost(.self, box = TRUE)  # If geomc is line, only box intersects make sense
    } else {
      boostGrid <- boost(.self)  # If geomc is line, only box intersects make sense
    }
    geomc.boost <- boost(geomc)
    intrs.ls <- bg_intersects(geomc.boost, boostGrid)

    # Iterate over polygons, bands, apply fun if not null
    missing.idx <- c()
    for (i in 1:length(intrs.ls)) {
      idx <- intrs.ls[[i]]
      if (length(idx) > 0) {
        if (is.null(fun)) {
          valmat <- matrix(NA, length(idx), nbands)
          for (band in 1:nbands) {
            valmat[,band] <- rasterbands[[band]][idx]
          }
          out[[i]] <- valmat
        } else {
          for (band in 1:nbands) {
            out[i,band][[1]] <- fun(rasterbands[[band]][idx])
          }
        }
      } else {
        missing.idx <- c(missing.idx, i)
      }
    }

    # If small activated: Fill missings using box intersect
    if (small & length(missing.idx) > 0 & !isLine) {

      # Create box grid, boost geometries, intersect
      boostBoxGrid <- boost(.self, box = TRUE)
      missing.boost <- geomc.boost[missing.idx]
      intrs.ls <- bg_intersects(missing.boost, boostBoxGrid)

      # Iterate over polygons, bands, apply fun if not null
      for (i in 1:length(intrs.ls)) {
        geom.idx <- missing.idx[i]
        idx <- intrs.ls[[i]]
        if (length(idx) > 0) {
          if (is.null(fun)) {
            valmat <- matrix(NA, length(idx), nbands)
            for (band in 1:nbands) {
              valmat[,band] <- rasterbands[[band]][idx]
            }
            out[[geom.idx]] <- valmat
          } else {
            for (band in 1:nbands) {
              out[geom.idx,band][[1]] <- fun(rasterbands[[band]][idx])
            }
          }
        }
      }
    }

  return(out)
})
