#' Generate inverse distance-weighted rasters
#'
#' \code{hydroweight::hydroweight()} generates distance-weighted rasters for targets
#' on a digital elevation model raster. Examples of targets include single points,
#' areas such as lakes, or linear features such as streams. The function outputs
#' a list of \code{length(weighting_scheme)} and an accompanying \code{*.rds} file
#' of distance-weighted rasters for targets (\code{target_O} is a point/area target as
#' in iFLO and \code{target_S} is a stream/waterbody target as in iFLS in
#' Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function acts on a single set of targets but can produce
#' multiple weights. The distance-weighted rasters, can be used for generating
#' distance-weighted using \code{hydroweight::hydroweight_attributes()} (e.g., % urban
#' cover weighted by flow distance to a point). See vignette for workflows.
#'
#' Spatial layers should align (i.e., identical coordinate reference systems - CRS).
#' Through processing, targets are converted to the resolution and CRS of the
#' digital elevation model/flow accumulation CRS. TEMP-* files are generated in
#'  \code{hydroweight_dir} depending on processing step and can be overwritten.
#'
#' For \code{weighting_scheme}:
#'
#' "lumped" indicates all weights = 1
#'
#' "iEucO" indicates Euclidean distance to \code{target_O}
#'
#' "iEucS" indicates Euclidean distance to \code{target_S}
#'
#' "iFLO" indicates d8 flow-path distance to \code{target_O}
#'
#' "iFLS" indicates d8 flow-path distance to \code{target_S}
#'
#' "HAiFLO" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to \code{target_O}
#'
#' "HAiFLS" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to \code{target_S}
#'
#' @param hydroweight_dir character. File path for read/write.
#' @param target_O \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "target.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Target for iEucO or iFLO.
#' @param target_S \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "target.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Target for iEucS or iFLS.
#' @param target_uid character. Unique identifier to precede exported list \code{*.rds} (i.e., "target_uid"_inv_distances.rds)
#' @param OS_combine logical. Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use \code{TRUE} or \code{FALSE}.
#' @param clip_region numeric, \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "clip_region.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Region over which distances are calculated. If numeric, \code{sf::sf_buffer()} produces a default buffer of crs-specific numeric width around \code{target_O}, exports, then clips using \code{whitebox}; if \code{sf}, exports and clips using \code{whitebox}; and if character, loads the file, converts to \code{sf}, exports, and clips using \code{whitebox}.
#' @param dem character (with extension, e.g., "dem.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Digital elevation model raster.
#' @param flow_accum  character (with extension, e.g., "flow_accum.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Flow accumulation raster (units: # of cells).
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param inv_function function. Inverse function used in \code{raster::calc()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1}
#' @return A named list of distance-weighted rasters and accompanying \code{*.rds} in \code{hydroweight_dir}
#' @export
#'
hydroweight <- function(hydroweight_dir = NULL,
                        target_O = NULL,
                        target_S = NULL,
                        target_uid = NULL,
                        OS_combine = NULL,
                        clip_region = NULL,
                        dem = NULL,
                        dem_crs = NULL,
                        flow_accum = NULL,
                        weighting_scheme = NULL,
                        inv_function = NULL) {
  message("Preparing hydroweight layers @ ", Sys.time())

  ## PREPARE HYDROWEIGHT LAYERS ----

  (dem_r <- raster::raster(file.path(hydroweight_dir, dem)))

  if(is.na(crs(dem_r)) | is.null(crs(dem_r))){
    stop("dem crs() is NULL or NA. Apply projection before continuing")
  }

  (dem_crs <- raster::crs(dem_r))

  ## Prepare clip_region ----
  if (!is.null(clip_region)) {
    if (class(clip_region)[1] == "numeric" & class(target_O)[1] == "sf") {
      clip_region <- sf::st_buffer(target_O, clip_region)
      sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
        append = FALSE, quiet = TRUE
      )
    }

    if (class(clip_region)[1] == "sf") {
      sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
        append = FALSE, quiet = TRUE
      )
    }

    if (class(clip_region)[1] == "character") {
      if (!grepl(".shp", clip_region) & !grepl(".tif", clip_region)) {
        stop("If clip_region is character, target_O should be a .shp or .tif")
      }

      if (grepl(".shp", clip_region)) {
        clip_region <- sf::st_read(file.path(hydroweight_dir, clip_region))
        sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
          append = FALSE, quiet = TRUE
        )
      }

      if (grepl(".tif", clip_region)) {
        clip_region <- raster::raster(file.path(hydroweight_dir, clip_region))

        if(is.na(crs(clip_region)) | is.null(crs(clip_region))){
          stop("clip_region crs() is NULL or NA. Apply projection before continuing")
        }

        whitebox::wbt_reclass(
          input = file.path(hydroweight_dir, clip_region),
          output = file.path(hydroweight_dir, "TEMP-clip_region.tif"),
          reclass_vals = "1;0.0,max")

        whitebox::wbt_raster_to_vector_polygons(
          input = file.path(hydroweight_dir, "TEMP-clip_region.tif"),
          output = file.path(hydroweight_dir, "TEMP-clip_region.shp")
        )

        #clip_region[!is.na(clip_region)] <- 1
        #clip_region <- raster::rasterToPolygons(clip_region, dissolve = TRUE)
        #clip_region <- sf::st_as_sf(clip_region)
        #sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
        #  append = FALSE, quiet = TRUE
        #)
      }
    }
  }

  if (is.null(clip_region)) {

    whitebox::wbt_reclass(
      input = file.path(hydroweight_dir, dem),
      output = file.path(hydroweight_dir, "TEMP-clip_region.tif"),
      reclass_vals = "1;0.0;30000", # because max sometimes misses max values?
      #verbose_mode = TRUE
      )

    whitebox::wbt_raster_to_vector_polygons(
      input = file.path(hydroweight_dir, "TEMP-clip_region.tif"),
      output = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
      #verbose_mode = TRUE
    )

    clip_region <- st_read(file.path(hydroweight_dir, "TEMP-clip_region.shp"), quiet = TRUE)
    st_crs(clip_region) <- st_crs(dem_crs)
    st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"), append = FALSE, quiet = TRUE)

    #clip_region <- raster::raster(file.path(hydroweight_dir, dem))
    #clip_region[!is.na(clip_region)] <- 1
    #clip_region <- raster::rasterToPolygons(clip_region, dissolve = TRUE)
    #clip_region <- sf::st_as_sf(clip_region)
    #sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
    #  append = FALSE, quiet = TRUE
    #)
  }

  whitebox::wbt_clip_raster_to_polygon(
    input = file.path(hydroweight_dir, dem),
    polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
    output = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
    verbose = TRUE
  )

  dem_clip <- raster::raster(file.path(hydroweight_dir, "TEMP-dem_clip.tif"))
  raster::crs(dem_clip) <- dem_crs

  ## Prepare target_O ----

  ## If sf, write to .shp
  if (class(target_O)[1] == "sf") {
    sf::st_write(target_O, file.path(hydroweight_dir, "TEMP-target_O.shp"),
      append = FALSE, quiet = TRUE
    )
  }

  ## If character, check if .shp or .tif, stop if not, then notify.
  ## If character and .shp or .tif, load and write to .shp
  if (class(target_O)[1] == "character") {
    if (!grepl(".shp", target_O) & !grepl(".tif", target_O)) {
      stop("If target_O is character, target_O should be a .shp or .tif")
    }

    if (grepl(".shp", target_O)) {
      target_O <- sf::st_read(file.path(hydroweight_dir, target_O))
      sf::st_write(target_O, file.path(hydroweight_dir, "TEMP-target_O.shp"),
        append = FALSE, quiet = TRUE
      )
    }

    if (grepl(".tif", target_O)) {
      target_O <- raster::raster(file.path(hydroweight_dir, target_O))
    }
  }

  ## write target_O and adjust to clip_region if sf
  if (class(target_O)[1] == "sf") {
    if (any(sf::st_is(target_O, c("POLYGON", "MULTIPOLYGON")))) {
      target_O_r <- raster::rasterize(sf::as_Spatial(target_O),
        y = dem_clip,
        filename = file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
        field = 1,
        overwrite = TRUE
      )
    }

    if (any(sf::st_is(target_O, c("LINESTRING", "MULTILINESTRING")))) {
      target_O_r <- raster::rasterize(sf::as_Spatial(target_O),
        y = dem_clip,
        filename = file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
        field = 1,
        overwrite = TRUE
      )
    }

    if (sf::st_is(target_O, "POINT")) {
      target_O_r <- raster::rasterize(sf::st_coordinates(target_O),
        y = dem_clip,
        field = 1,
        filename = file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
        overwrite = TRUE
      )
    }
  }

  ## write target_O and adjust to clip_region if RasterLayer
  if (class(target_O)[1] == "RasterLayer") {
    target_O_r <- target_O
    target_O_r <- raster::projectRaster(target_O_r, dem_clip, method = "ngb")
    raster::writeRaster(target_O_r, file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
      overwrite = TRUE
    )
  }

  ## Prepare target_S ----

  ## If sf, write to .shp

  if (class(target_S)[1] == "sf") {
    sf::st_write(target_S, file.path(hydroweight_dir, "TEMP-target_S.shp"),
      append = FALSE, quiet = TRUE
    )
  }

  ## If character, check if .shp or .tif, stop if not, then notify.
  ## If character and .shp or .tif, load and write to .shp
  if (class(target_S)[1] == "character") {
    if (!grepl(".shp", target_S) & !grepl(".tif", target_S)) {
      stop("If target_S is character, target_S should be a .shp or .tif")
    }

    if (grepl(".shp", target_S)) {
      target_S <- sf::st_read(file.path(hydroweight_dir, target_S))
      sf::st_write(target_S, file.path(hydroweight_dir, "TEMP-target_S.shp"),
        append = FALSE, quiet = TRUE
      )
    }

    if (grepl(".tif", target_S)) {
      target_S <- raster::raster(file.path(hydroweight_dir, target_S))
    }
  }

  ## write target_S and adjust to clip_region if sf
  if (class(target_S)[1] == "sf") {

    if (any(sf::st_is(target_S, c("POLYGON", "MULTIPOLYGON")))) {
      target_S_r <- raster::rasterize(sf::as_Spatial(target_S),
        y = dem_clip,
        filename = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        field = 1,
        overwrite = TRUE
      )
    }

    if (any(sf::st_is(target_S, c("LINESTRING", "MULTILINESTRING")))) {
      target_S_r <- raster::rasterize(sf::as_Spatial(target_S),
        y = dem_clip,
        filename = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        field = 1,
        overwrite = TRUE
      )
    }

    if (sf::st_is(target_S, "POINT")) {
      target_S_r <- raster::rasterize(sf::st_coordinates(target_S),
        y = dem_clip,
        field = 1,
        filename = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        overwrite = TRUE
      )
    }
  }

  ## write target_S and adjust to clip_region if RasterLayer
  if (class(target_S)[1] == "RasterLayer") {
    target_S_r <- target_S
    target_S_r <- raster::projectRaster(target_S_r, dem_clip, method = "ngb")
    raster::writeRaster(target_S_r, file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
      overwrite = TRUE
    )
  }

  ## Prepare OS_combine ----
  if (is.null(OS_combine)) {
    OS_combine <- FALSE
  }

  if (OS_combine == TRUE){

  (target_S_OS <- target_S_r)
  (target_O_OS <- target_O_r)

  target_S_OS[target_S_OS > 0] <- 1
  target_O_OS[target_O_OS > 0] <- 1

  mosaic_list <- list(target_O_OS, target_S_OS)
  mosaic_list$fun <- sum
  mosaic_list$na.rm <- TRUE

  OS_combine_r <- do.call(raster::mosaic, mosaic_list)
  OS_combine_r[OS_combine_r > 0] <- 1
  OS_combine_r[OS_combine_r == 0] <- -9999

  raster::crs(OS_combine_r) <- dem_crs
  raster::writeRaster(OS_combine_r,
                      file.path(hydroweight_dir, "TEMP-OS_combine.tif"),
                      overwrite = TRUE, options = c("COMPRESS=NONE"),
                      Naflag = -9999
  )
  }

  ## RUN DISTANCE-WEIGHTING ----

  message("Running distance-weighting @ ", Sys.time())

  if ("lumped" %in% weighting_scheme) {

    whitebox::wbt_reclass(
      input = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
      output = file.path(hydroweight_dir, "TEMP-lumped.tif"),
<<<<<<< HEAD
      reclass_vals = "1;0.0;30000")  # because max sometimes misses max values?

    lumped_inv <- raster::raster(file.path(hydroweight_dir, "TEMP-lumped.tif"), values = TRUE)
    raster::crs(lumped_inv) <- dem_crs
    lumped_inv <- setValues(raster(lumped_inv), lumped_inv[])
=======
      reclass_vals = "1;0.0;max")

    lumped_inv <- raster::raster(file.path(hydroweight_dir, "TEMP-lumped.tif"))
    raster::crs(lumped_inv) <- dem_crs

>>>>>>> 911612189952d5d2b74019afb6b2922319d311d2
    #lumped_inv <- dem_clip
    #lumped_inv[!is.na(lumped_inv)] <- 1
  }

  ## iEucO, Euclidean distance to target_O ----
  if ("iEucO" %in% weighting_scheme) {

    whitebox::wbt_reclass(
      input = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
      output = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
<<<<<<< HEAD
      reclass_vals = "1;0.0;30000")  # because max sometimes misses max values?
=======
      reclass_vals = "1;0.0;max")
>>>>>>> 911612189952d5d2b74019afb6b2922319d311d2

    #cost <- dem_clip
    #cost[cost > 0] <- 1
    #raster::writeRaster(cost,
    #  file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
    #  overwrite = TRUE, options = c("COMPRESS=NONE")
    #)

    whitebox::wbt_cost_distance(
      source = file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
      cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
      out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
      out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
      verbose_mode = TRUE
    )

    iEucO <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
    iEucO_inv <- raster::calc(iEucO, fun = inv_function)
    raster::crs(iEucO_inv) <- dem_crs

    raster::writeRaster(iEucO_inv,
      file.path(hydroweight_dir, "TEMP-iEucO.tif"),
      overwrite = TRUE, options = c("COMPRESS=NONE"),
      Naflag = -9999
    )
  }

  ## iEucS, Euclidean distance to streams ----
  if ("iEucS" %in% weighting_scheme) {

    if (OS_combine == FALSE) {

      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
        out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
        out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
        verbose_mode = TRUE
      )

      iEucS <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
      iEucS_inv <- raster::calc(iEucS, fun = inv_function)
      raster::crs(iEucS_inv) <- dem_crs

      raster::writeRaster(iEucS_inv,
        file.path(hydroweight_dir, "TEMP-iEucS.tif"),
        overwrite = TRUE, options = c("COMPRESS=NONE"),
        Naflag = -9999
      )
    }

    if (OS_combine == TRUE) {

      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, "TEMP-OS_combine.tif"),
        cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
        out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
        out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
        verbose_mode = TRUE
      )

      iEucS <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
      iEucS_inv <- raster::calc(iEucS, fun = inv_function)
      raster::crs(iEucS_inv) <- dem_crs

      raster::writeRaster(iEucS_inv,
        file.path(hydroweight_dir, "TEMP-iEucS.tif"),
        overwrite = TRUE, options = c("COMPRESS=NONE"),
        Naflag = -9999
      )
    }
  }

  ## iFLO, flow line distance to target_O ----
  if ("iFLO" %in% weighting_scheme) {

    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
      streams = file.path(hydroweight_dir, "TEMP-target_O_clip.tif"),
      output = file.path(hydroweight_dir, "TEMP-flowdist.tif"),
      verbose_mode = TRUE
    )

    iFLO <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
    iFLO_inv <- raster::calc(iFLO, fun = inv_function)
    raster::crs(iFLO_inv) <- dem_crs

    raster::writeRaster(iFLO_inv,
      file.path(hydroweight_dir, "TEMP-iFLO.tif"),
      overwrite = TRUE, options = c("COMPRESS=NONE"),
      Naflag = -9999
    )
  }

  ## iFLS, flow line distance to target_S ----
  if ("iFLS" %in% weighting_scheme) {
    if (OS_combine == FALSE) {

      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
        streams = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        output = file.path(hydroweight_dir, "TEMP-flowdist.tif"),
        verbose_mode = TRUE
      )

      iFLS <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
      iFLS_inv <- raster::calc(iFLS, fun = inv_function)
      raster::crs(iFLS_inv) <- dem_crs

      raster::writeRaster(iFLS_inv,
        file.path(hydroweight_dir, "TEMP-iFLS.tif"),
        overwrite = TRUE, options = c("COMPRESS=NONE"),
        Naflag = -9999
      )
    }

    if (OS_combine == TRUE) {

      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
        streams = file.path(hydroweight_dir, "TEMP-OS_combine.tif"),
        output = file.path(hydroweight_dir, "TEMP-flowdist.tif"),
        verbose_mode = TRUE
      )

      iFLS <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
      iFLS_inv <- raster::calc(iFLS, fun = inv_function)
      raster::crs(iFLS_inv) <- dem_crs

      raster::writeRaster(iFLS_inv,
        file.path(hydroweight_dir, "TEMP-iFLS.tif"),
        overwrite = TRUE, options = c("COMPRESS=NONE"),
        Naflag = -9999
      )
    }
  }

  ## HA-iFLO and HA-iFLS, hydrologically active ----
  if (("HAiFLO" %in% weighting_scheme) | ("HAiFLS" %in% weighting_scheme)) {

    whitebox::wbt_clip_raster_to_polygon(
      input = file.path(hydroweight_dir, flow_accum),
      polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
      output = file.path(hydroweight_dir, "TEMP-flow_accum_clip.tif"),
      verbose = TRUE
    )

    accum_clip <- raster::raster(file.path(hydroweight_dir, "TEMP-flow_accum_clip.tif"))
    accum_clip <- accum_clip + 1
    raster::crs(accum_clip) <- dem_crs

    if (OS_combine == TRUE) {
      if ("HAiFLO" %in% weighting_scheme) {

        HAiFLO_inv <- iFLO_inv * accum_clip
        HAiFLO_inv <- raster::mask(HAiFLO_inv, OS_combine_r, maskvalue = 1)

        raster::writeRaster(HAiFLO_inv,
          file.path(hydroweight_dir, "TEMP-HAiFLO.tif"),
          overwrite = TRUE, options = c("COMPRESS=NONE"),
          Naflag = -9999
        )
      }

      if ("HAiFLS" %in% weighting_scheme) {

        HAiFLS_inv <- iFLS_inv * accum_clip
        HAiFLS_inv <- raster::mask(HAiFLS_inv, OS_combine_r, maskvalue = 1)

        raster::writeRaster(HAiFLS_inv,
          file.path(hydroweight_dir, "TEMP-HAiFLS.tif"),
          overwrite = TRUE, options = c("COMPRESS=NONE"),
          Naflag = -9999
        )
      }
    }

    if (OS_combine == FALSE) {
      if ("HAiFLO" %in% weighting_scheme) {

        HAiFLO_inv <- iFLO_inv * accum_clip
        HAiFLO_inv <- raster::mask(HAiFLO_inv, target_O_r, maskvalue = 1)

        raster::writeRaster(HAiFLO_inv,
          file.path(hydroweight_dir, "TEMP-HAiFLO.tif"),
          overwrite = TRUE, options = c("COMPRESS=NONE"),
          Naflag = -9999
        )
      }

      if ("HAiFLS" %in% weighting_scheme) {

        HAiFLS_inv <- iFLS_inv * accum_clip
        HAiFLS_inv <- raster::mask(HAiFLS_inv, target_S_r, maskvalue = 1)

        raster::writeRaster(HAiFLS_inv,
          file.path(hydroweight_dir, "TEMP-HAiFLS.tif"),
          overwrite = TRUE, options = c("COMPRESS=NONE"),
          Naflag = -9999
        )
      }
    }
  }


  ## PREPARE OUTPUT ----

  weighting_scheme_inv <- paste0(weighting_scheme, "_inv")

  dist_list <- lapply(weighting_scheme_inv, function(x) {
    get(x)
  })
  names(dist_list) <- weighting_scheme

  saveRDS(dist_list, file = file.path(
    hydroweight_dir,
    paste0(target_uid, "_inv_distances.rds")
  ))

  return(dist_list)
}
