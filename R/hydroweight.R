#' Generate inverse distance-weighted rasters
#'
#' \code{hydroweight::hydroweight()} generates distance-weighted rasters for targets
#' on a digital elevation model raster. Examples of targets include single points,
#' areas such as lakes, or linear features such as streams. The function outputs
#' a list of \code{length(weighting_scheme)} and an accompanying \code{*.rds} file
#' of distance-weighted rasters for targets (\code{target_O} is a point target as
#' in iFLO and \code{target_S} is a stream/waterbody target as in iFLS in
#' Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function produces one instance of the weighting (i.e., one target and multiple weights).
#' Distance-weighted rasters, in turn, can be used in producing distance-weighted
#' attributes using \code{hydroweight::hydroweight_attributes()} (e.g., % urban
#' cover weighted by flow distance to a point). See examples and vignette for simple workflow.
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
#' "EucO" indicates Euclidean distance to \code{target_O}
#'
#' "EucS" indicates Euclidean distance to \code{target_S}
#'
#' "FLO" indicates d8 flow-path distance to \code{target_O}
#'
#' "FLS" indicates d8 flow-path distance to \code{target_S}
#'
#' "HAFLO" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to \code{target_O}
#'
#' "HAFLS" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to \code{target_S}
#'
#' @param hydroweight_dir character. File path for read/write.
#' @param target_O \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "target.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Target for EucO or FLO.
#' @param target_S \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "target.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Target for EucS or FLS.
#' @param OS_combine logical. Should target_O and target_S be merged as targets for EucS, FLS, and/or HAFLS? Use \code{TRUE} or \code{FALSE}.
#' @param clip_region numeric, \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "clip_region.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Region over which distances are calculated. If numeric, \code{sf::sf_buffer()} produces a default buffer of crs-specific numeric width around \code{target_O}, exports, then clips using \code{whitebox}; if \code{sf}, exports and clips using \code{whitebox}; and if character, loads the file, converts to \code{sf}, exports, and clips using \code{whitebox}.
#' @param dem character (with extension, e.g., "dem.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Digital elevation model raster.
#' @param dem_crs \code{CRS-class}. Digital elevation model coordinate reference system.
#' @param flow_accum  character (with extension, e.g., "flow_accum.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Flow accumulation raster (units: # of cells).
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "EucO", "EucS", "FLO", "FLS", "HAFLO", "HAFLS")
#' @param inv_function function. Inverse function used in \code{raster::calc()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1}
#' @param uid character. Unique identifier to precede exported list \code{*.rds} (i.e., uid_hydrology.rds)
#' @return A named list of distance-weighted rasters and accompanying \code{*.rds} in \code{hydroweight_dir}
#' @export
#'
hydroweight <- function(hydroweight_dir = NULL,
                        target_O = NULL,
                        target_S = NULL,
                        OS_combine = NULL,
                        clip_region = NULL,
                        dem = NULL,
                        dem_crs = NULL,
                        flow_accum = NULL,
                        weighting_scheme = NULL,
                        inv_function = NULL,
                        uid = NULL
                        ){

  message("Preparing hydroweight layers @ ", Sys.time())

  ## PREPARE HYDROWEIGHT LAYERS ----

  ## Get dem crs if not specified
  if(is.null(dem_crs)){
    dem_crs <- raster::crs(
      raster::raster(file.path(hydroweight_dir, dem))
    )
  }

  ## write target_O to .shp ----

  ## If sf, write to .shp
  if (class(target_O)[1] == "sf") {
    sf::st_write(target_O, file.path(hydroweight_dir, "TEMP-target_O.shp"),
                 append = FALSE, quiet = TRUE
    )
  }

  ## If character, check if .shp or .tif, stop if not, then notify.
  ## If character and .shp or .tif, load and write to .shp
  if (class(target_O)[1] == "character") {

    if( !grepl(".shp", target_O) | !grepl(".tif", target_O) ){
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

  ## write clip_region to .shp ----

  if (class(clip_region)[1] == "numeric" & class(target_O)[1] == "sf") {
    clip_region <- sf::st_buffer(target_O, clip_region)
    sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
                 append = FALSE, quiet = TRUE
    )
  }

  if (class(clip_region)[1] == "sf") {
    sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
                 append = FALSE, quiet = T
    )
  }

  if (class(clip_region)[1] == "character") {

    if( !grepl(".shp", clip_region) | !grepl(".tif", clip_region) ){
      stop("If clip_region is character, target_O should be a .shp or .tif")
    }

    if (grepl(".shp", clip_region)) {
      clip_region <- sf::st_read(file.path(hydroweight_dir, clip_region))
      sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
                   append = FALSE
      )
    }

    if (grepl(".tif", clip_region)) {
      clip_region <- raster::raster(file.path(hydroweight_dir, clip_region))
      clip_region <- raster::rasterToPolygons(clip_region, dissolve = TRUE)
      clip_region <- sf::st_as_sf(clip_region)
      sf::st_write(clip_region, file.path(hydroweight_dir, "TEMP-clip_region.shp"),
                   append = FALSE
      )
    }
  }

  ## clip dem and target_S by clip_region ----

  whitebox::wbt_clip_raster_to_polygon(
    input = file.path(hydroweight_dir, dem),
    polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
    output = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
    verbose = TRUE
  )

  dem_clip <- raster::raster(file.path(hydroweight_dir, "TEMP-dem_clip.tif"))
  raster::crs(dem_clip) <- dem_crs

  ## write target_O to .tif ----
  if(class(target_O)[1] == "RasterLayer"){

    target_O_r <- target_O
    raster::crs(target_O_r) <- dem_crs
    raster::writeRaster(target_O_r, file.path(hydroweight_dir, "TEMP-target_O.tif"),
                        overwrite = TRUE
    )

  }

  if(class(target_O)[1] == "sf"){

    target_O_r <- fasterize::fasterize(target_O, raster = dem_clip)
    raster::crs(target_O_r) <- dem_crs
    raster::writeRaster(target_O_r, file.path(hydroweight_dir, "TEMP-target_O.tif"),
                        overwrite = TRUE
    )

  }

  ## write target_S to .tif ----
  if (class(target_S)[1] == "sf") {
    sf::st_write(target_S, file.path(hydroweight_dir, "TEMP-target_S.shp"),
                 append = FALSE, quiet = TRUE
    )

    ## Convert target to raster
    target_S_r <- fasterize::fasterize(target_S, raster = dem_clip)
    raster::crs(target_S_r) <- dem_crs
    raster::writeRaster(target_S_r, file.path(hydroweight_dir, "TEMP-target_S.tif"),
                        overwrite = TRUE
    )
  }

  if(class(target_S)[1] == "RasterLayer"){

    target_S_r <- target_S
    raster::crs(target_S_r) <- dem_crs
    raster::writeRaster(target_S_r, file.path(hydroweight_dir, "TEMP-target_S.tif"),
                        overwrite = TRUE
    )

    whitebox::wbt_clip_raster_to_polygon(
      input = file.path(hydroweight_dir, "TEMP-target_S.tif"),
      polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
      output = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
      verbose = TRUE
    )

  }

  if (class(target_S)[1] == "character") {

    if( !grepl(".shp", target_S) | !grepl(".tif", target_S) ){
      stop("If target_S is character, target_S should be a .shp or .tif")
    }

    if (grepl(".shp", target_S)) {
      target_S <- sf::st_read(file.path(hydroweight_dir, target_S))
      sf::st_write(target_S, file.path(hydroweight_dir, "TEMP-target_S.shp"),
                   append = FALSE
      )

      ## Convert target to raster
      target_S <- fasterize::fasterize(target_S, raster = dem_clip)
      raster::crs(target_S) <- dem_crs
      raster::writeRaster(target_S, file.path(hydroweight_dir, "TEMP-target_S.tif"),
                          overwrite = TRUE
      )
    }

    if (grepl(".tif", target_S)) {

      whitebox::wbt_clip_raster_to_polygon(
        input = file.path(hydroweight_dir, target_S),
        polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
        output = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        verbose = TRUE
      )

    }
  }

  target_S_clip <- raster::raster(file.path(hydroweight_dir, "TEMP-target_S_clip.tif"))
  raster::crs(target_S_clip) <- dem_crs

  ## RUN DISTANCE-WEIGHTING ----

  message("Running distance-weighting @ ", Sys.time())

  ## allows user to ignore OS_combine if only doing stream distances
  if (is.null(OS_combine)) { OS_combine <- FALSE }

  if ("lumped" %in% weighting_scheme){

    lumped_inv <- dem_clip
    lumped_inv[!is.na(lumped_inv)] <- 1
  }

  ## EucO, Euclidean distance to target_O ----
  if ("EucO" %in% weighting_scheme) {

    cost <- dem_clip
    cost[cost > 0] <- 1
    raster::writeRaster(cost,
                        file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
                        overwrite = TRUE, options = c("COMPRESS=NONE")
    )

    whitebox::wbt_cost_distance(
      source = file.path(hydroweight_dir, "TEMP-target_O.tif"),
      cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
      out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
      out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
      verbose_mode = TRUE
    )

    EucO <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
    EucO_inv <- raster::calc(EucO, fun = inv_function)

    raster::writeRaster(EucO_inv,
                        file.path(hydroweight_dir, "TEMP-EucO.tif"),
                        overwrite = TRUE, options = c("COMPRESS=NONE"),
                        Naflag = -9999
    )

  }

  ## EucS, Euclidean distance to streams ----
  if ("EucS" %in% weighting_scheme) {

    if (OS_combine == FALSE) {
      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
        out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
        out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
        verbose_mode = TRUE
      )

      EucS <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
      EucS_inv <- raster::calc(EucS, fun = inv_function)

      raster::writeRaster(EucS_inv,
                          file.path(hydroweight_dir, "TEMP-EucS.tif"),
                          overwrite = TRUE, options = c("COMPRESS=NONE"),
                          Naflag = -9999
      )
    }

    if (OS_combine == TRUE) {
      target_S_OS <- target_S_clip
      target_O_OS <- target_O_r

      target_S_OS[target_S_OS > 0] <- 1
      target_O_OS[target_O_OS > 0] <- 1

      mosaic_list <- list(target_O_OS, target_S_OS)
      mosaic_list$fun <- sum
      mosaic_list$na.rm <- TRUE

      OS_combine_r <- do.call(mosaic, mosaic_list)
      OS_combine_r[OS_combine_r > 0] <- 1
      OS_combine_r[OS_combine_r == 0] <- -9999

      raster::crs(OS_combine_r) <- dem_crs
      raster::writeRaster(OS_combine_r,
                          file.path(hydroweight_dir, "TEMP-OS_combine.tif"),
                          overwrite = TRUE, options = c("COMPRESS=NONE"),
                          Naflag = -9999
      )

      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, "TEMP-OS_combine.tif"),
        cost = file.path(hydroweight_dir, "TEMP_dem_clip_cost.tif"),
        out_accum = file.path(hydroweight_dir, "TEMP-cost_distance.tif"),
        out_backlink = file.path(hydroweight_dir, "TEMP-cost_backlink.tif"),
        verbose_mode = TRUE
      )

      EucS <- raster::raster(file.path(hydroweight_dir, "TEMP-cost_distance.tif"))
      EucS_inv <- raster::calc(EucS, fun = inv_function)

      raster::writeRaster(EucS_inv,
                          file.path(hydroweight_dir, "TEMP-EucS.tif"),
                          overwrite = TRUE, options = c("COMPRESS=NONE"),
                          Naflag = -9999
      )
    }

  }

  ## FLO, flow line distance to target_O ----
  if ("FLO" %in% weighting_scheme) {

    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
      streams = file.path(hydroweight_dir, "TEMP-target_O.tif"),
      output = file.path(hydroweight_dir, "TEMP-flowdist.tif"),
      verbose_mode = TRUE
    )

    FLO <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
    FLO_inv <- raster::calc(FLO, fun = inv_function)
    raster::crs(FLO_inv) <- dem_crs

    raster::writeRaster(FLO_inv,
                        file.path(hydroweight_dir, "TEMP-FLO.tif"),
                        overwrite = TRUE, options = c("COMPRESS=NONE"),
                        Naflag = -9999
    )

  }

  ## FLS, flow line distance to target_S ----
  if ("FLS" %in% weighting_scheme) {

    if (OS_combine == FALSE) {
      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(hydroweight_dir, "TEMP-dem_clip.tif"),
        streams = file.path(hydroweight_dir, "TEMP-target_S_clip.tif"),
        output = file.path(hydroweight_dir, "TEMP-flowdist.tif"),
        verbose_mode = TRUE
      )

      FLS <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
      FLS_inv <- raster::calc(FLS, fun = inv_function)
      raster::crs(FLS_inv) <- dem_crs

      raster::writeRaster(FLS_inv,
                          file.path(hydroweight_dir, "TEMP-FLS.tif"),
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

      FLS <- raster::raster(file.path(hydroweight_dir, "TEMP-flowdist.tif"))
      FLS_inv <- raster::calc(FLS, fun = inv_function)
      raster::crs(FLS_inv) <- dem_crs

      raster::writeRaster(FLS_inv,
                          file.path(hydroweight_dir, "TEMP-FLS.tif"),
                          overwrite = TRUE, options = c("COMPRESS=NONE"),
                          Naflag = -9999
      )
    }
  }

  ## HA-FLO and HA-FLS, hydrologically active ----
  if (("HAFLO" %in% weighting_scheme) | ("HAFLS" %in% weighting_scheme)) {

    whitebox::wbt_clip_raster_to_polygon(
      input = file.path(hydroweight_dir, flow_accum),
      polygons = file.path(hydroweight_dir, "TEMP-clip_region.shp"),
      output = file.path(hydroweight_dir, "TEMP-flow_accum_clip.tif"),
      verbose = TRUE
    )

    accum_clip <- raster::raster(file.path(hydroweight_dir, "TEMP-flow_accum_clip.tif"))

    if (OS_combine == TRUE) {

      if ("HAFLO" %in% weighting_scheme) {
        HAFLO_inv <- FLO_inv * accum_clip
        HAFLO_inv <- raster::mask(HAFLO_inv, OS_combine_r, maskvalue = 1)

        raster::writeRaster(HAFLO_inv,
                            file.path(hydroweight_dir, "TEMP-HAFLO.tif"),
                            overwrite = TRUE, options = c("COMPRESS=NONE"),
                            Naflag = -9999
        )
      }

      if ("HAFLS" %in% weighting_scheme) {
        HAFLS_inv <- FLS_inv * accum_clip
        HAFLS_inv <- raster::mask(HAFLS_inv, OS_combine_r, maskvalue = 1)

        raster::writeRaster(HAFLS_inv,
                            file.path(hydroweight_dir, "TEMP-HAFLS.tif"),
                            overwrite = TRUE, options = c("COMPRESS=NONE"),
                            Naflag = -9999
        )
      }
    }

    if (OS_combine == FALSE) {

      if ("HAFLO" %in% weighting_scheme) {
        HAFLO_inv <- FLO_inv * accum_clip
        HAFLO_inv <- raster::mask(HAFLO_inv, target_S_clip, maskvalue = 1)

        raster::writeRaster(HAFLO_inv,
                            file.path(hydroweight_dir, "TEMP-HAFLO.tif"),
                            overwrite = TRUE, options = c("COMPRESS=NONE"),
                            Naflag = -9999
        )

      }

      if ("HAFLS" %in% weighting_scheme) {
        HAFLS_inv <- FLS_inv * accum_clip
        HAFLS_inv <- raster::mask(HAFLS_inv, target_S_clip, maskvalue = 1)

        raster::writeRaster(HAFLS_inv,
                            file.path(hydroweight_dir, "TEMP-HAFLS.tif"),
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
  names(dist_list) <- weighting_scheme_inv

  saveRDS(dist_list, file = file.path(
    hydroweight_dir,
    paste0(uid, "_inv_distances.rds")
  ))

  return(dist_list)
}
