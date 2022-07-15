#' Generate inverse distance-weighted rasters
#'
#' \code{hydroweight::hydroweight()} generates distance-weighted rasters for targets
#' on a digital elevation model raster. Examples of targets include single points,
#' areas such as lakes, or linear features such as streams. The function outputs
#' a list of \code{length(weighting_scheme)} and an accompanying \code{*.rds} file
#' of distance-weighted rasters for targets (\code{target_O} is a point/area target as
#' in iFLO and \code{target_S} is a stream/linear feature target as in iFLS in
#' Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function acts on a single set of targets but can produce
#' multiple weights. The distance-weighted rasters, can be used for generating
#' distance-weighted landscape statistics using \code{hydroweight::hydroweight_attributes()} (e.g., % urban
#' cover weighted by flow distance to a point). See https://github.com/bkielstr/hydroweight for workflows.
#'
#' Spatial layers should align (i.e., identical coordinate reference systems - CRS).
#' Through processing, targets are converted to the resolution and CRS of the
#' digital elevation model/flow accumulation CRS. _TEMP-* files are generated in
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
#' @param OS_combine logical. Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use \code{TRUE} or \code{FALSE}. This allows cells surrounding \code{target_O} to flow directly into \code{target_O} rather than be forced through \code{target_S}.
#' @param clip_region numeric, \code{sf}, \code{RasterLayer}, or character (with extension, e.g., "clip_region.shp") of file found in \code{hydroweight_dir} of ESRI Shapefile type or GeoTiFF type only. Region over which distances are calculated. If numeric, \code{sf::sf_buffer()} produces a default buffer of crs-specific numeric width around \code{target_O}, exports, then clips using \code{whitebox}; if \code{sf}, exports and clips using \code{whitebox}; and if character, loads the file, converts to \code{sf}, exports, and clips using \code{whitebox}.
#' @param dem character (with extension, e.g., "dem.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Digital elevation model raster.
#' @param flow_accum  character (with extension, e.g., "flow_accum.tif") of file found in \code{hydroweight_dir} of GeoTiFF type. Flow accumulation raster (units: # of cells).
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @return Named list of distance-weighted rasters and accompanying \code{*.rds} in \code{hydroweight_dir}
#' @export
#'
hydroweight <- function(hydroweight_dir = NULL,
                        target_O = NULL,
                        target_S = NULL,
                        target_uid = NULL,
                        OS_combine = NULL,
                        clip_region = NULL,
                        dem = NULL,
                        flow_accum = NULL,
                        weighting_scheme = NULL,
                        inv_function = function(x) {
                          (x * 0.001 + 1)^-1
                        }) {
  # TODO: Add checks to make sure all raster inputs are SpatRaster NOT RasterLayer
  # TODO: class(___)[1] should be changed to inherits()

  if (is.null(target_uid)) stop("target_uid must be specified")
  ## SET UP ----

  ## Set up raster temp file for out-of-memory raster files
  terra::terraOptions(tempdir = file.path(tempdir(),basename(tempfile())))
  #terra::rastTmpFile(prefix = paste0(target_uid,"_hydroweight_")) #PS: this isn't an option in terra

  ## Set whitebox verbose_mode to FALSE
  whitebox::wbt_options(verbose = FALSE)

  ## PREPARE HYDROWEIGHT LAYERS ----

  message("Preparing hydroweight layers @ ", Sys.time())

  (dem_r <- terra::rast(file.path(hydroweight_dir, dem)))

  if (is.na(terra::crs(dem_r)) | is.null(terra::crs(dem_r))) {
    stop("dem crs() is NULL or NA. Apply projection before continuing")
  }

  (dem_crs <- terra::crs(dem_r))

  ## Prepare clip_region ----
  clip_region_path<-file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.shp"))
  if (!is.null(clip_region)) {
    if (class(clip_region)[1] == "numeric" & class(target_O)[1] == "sf") {
      clip_region <- sf::st_buffer(target_O, clip_region)
      sf::st_write(clip_region, clip_region_path,append = FALSE, quiet = TRUE
      )
    }

    if (class(clip_region)[1] == "sf") {
      sf::st_write(clip_region, clip_region_path, append = FALSE, quiet = TRUE
      )
    }

    if (class(clip_region)[1] == "character") {
      if (!grepl(".shp", clip_region) & !grepl(".tif", clip_region)) {
        stop("If clip_region is character, target_O should be a .shp or .tif")
      }

      if (grepl(".shp", clip_region)) {
        clip_region <- sf::st_read(file.path(hydroweight_dir, clip_region))
        sf::st_write(clip_region, clip_region_path,append = FALSE, quiet = TRUE
        )
      } else if (grepl(".tif", clip_region)) {
        #clip_region <- terra::rast(file.path(hydroweight_dir, clip_region))

        #if (is.na(terra::crs(clip_region)) | is.null(terra::crs(clip_region))) {
        #  stop("clip_region crs() is NULL or NA. Apply projection before continuing")
        #}

        whitebox::wbt_reclass(
          input = file.path(hydroweight_dir, clip_region),
          output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.tif")),
          reclass_vals = "1,min,max"
        )

        whitebox::wbt_raster_to_vector_polygons(
          input = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.tif")),
          output = clip_region_path
        )

        clip_region <- sf::st_read(clip_region_path, quiet = TRUE)
        sf::st_crs(clip_region) <- sf::st_crs(dem_crs)
        sf::st_write(clip_region, clip_region_path, append = FALSE, quiet = TRUE)
      }
    }


  }

  if (is.null(clip_region)) {

    whitebox::wbt_reclass(
      input = file.path(hydroweight_dir, dem),
      output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.tif")),
      reclass_vals = "1,min,1000000", # because max sometimes misses max values?
      verbose_mode = FALSE
    )

    whitebox::wbt_raster_to_vector_polygons(
      input = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.tif")),
      output = clip_region_path,
      verbose_mode = FALSE
    )

    clip_region <- sf::st_read(clip_region_path, quiet = TRUE)
    sf::st_crs(clip_region) <- sf::st_crs(dem_crs)
    sf::st_write(clip_region, clip_region_path, append = FALSE, quiet = TRUE)

    # clip_region <- terra::rast(file.path(hydroweight_dir, dem))
    # clip_region[!is.na(clip_region)] <- 1
    # clip_region <- terra::rastToPolygons(clip_region, dissolve = TRUE)
    # clip_region <- sf::st_as_sf(clip_region)
    # sf::st_write(clip_region, file.path(hydroweight_dir, "_TEMP-clip_region.shp"),
    #  append = FALSE, quiet = TRUE
    # )
  }

  dem_clip_path<-file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif"))


  dem_clip<-terra::crop(
    x=terra::rast(file.path(hydroweight_dir, dem)),
    y=terra::vect(clip_region_path)
  )

  dem_clip<-terra::mask(
    x=dem_clip,
    mask=terra::vect(clip_region_path),
    touches=T
  )

  terra::crs(dem_clip)<-dem_crs
  terra::writeRaster(dem_clip,
                     dem_clip_path,
                     overwrite=T)
  # whitebox::wbt_clip_raster_to_polygon( #PS I think is is very slow for large rasters
  #   input = file.path(hydroweight_dir, dem),
  #   polygons = clip_region_path,
  #   output = dem_clip_path,
  #   verbose_mode = FALSE
  # )


  dem_clip <- terra::rast(dem_clip_path)
  terra::crs(dem_clip) <- dem_crs

  ## Prepare target_O ----

  ## If sf, write to .shp
  if (class(target_O)[1] == "sf") {
    sf::st_write(target_O, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O.shp")),
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
      sf::st_write(target_O, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O.shp")),
                   append = FALSE, quiet = TRUE
      )
    } else if (grepl(".tif", target_O)) {
      target_O <- terra::rast(file.path(hydroweight_dir, target_O))
    }
  }

  ## write target_O and adjust to clip_region if sf
  if (class(target_O)[1] == "sf") {
    if (any(sf::st_is(target_O, c("POLYGON", "MULTIPOLYGON")))) {
      target_O_r <- terra::rasterize(sf::as_Spatial(target_O),
                                     y = dem_clip,
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
                                     field = 1,
                                     overwrite = TRUE
      )
    }

    if (any(sf::st_is(target_O, c("LINESTRING", "MULTILINESTRING")))) {
      dem_clip_extent <- sf::st_as_sfc(sf::st_bbox(terra::ext(dem_clip)))
      sf::st_crs(dem_clip_extent) <- sf::st_crs(dem_clip)

      target_O_int <- sf::st_intersects(target_O, dem_clip_extent)
      target_O_int <- target_O[lengths(target_O_int) > 0, ]

      target_O_r <- terra::rasterize(sf::as_Spatial(target_O_int),
                                     y = dem_clip,
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
                                     field = 1,
                                     overwrite = TRUE
      )
    }

    if (sf::st_is(target_O, "POINT")) {
      target_O_r <- terra::rasterize(sf::st_coordinates(target_O),
                                     y = dem_clip,
                                     # field = 1, # this isn't needed in terra
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
                                     overwrite = TRUE
      )
    }
  }

  ## write target_O and adjust to clip_region if RasterLayer
  if (class(target_O)[1] %in% c("SpatRaster","RasterLayer")) {
    if (class(target_O)[1] %in% c("RasterLayer")) target_O_r <- terra::rast(target_O)
    if (class(target_O)[1] %in% c("SpatRaster")) target_O_r <- target_O

    target_O_r <- terra::project(target_O_r, dem_clip, method = "near",overwrite = TRUE,filename=file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")))
    terra::writeRaster(target_O_r, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
                       overwrite = TRUE
    )
  }

  ## Prepare target_S ----

  ## If sf, write to .shp

  if (class(target_S)[1] == "sf") {
    sf::st_write(target_S, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S.shp")),
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
      sf::st_write(target_S, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S.shp")),
                   append = FALSE, quiet = TRUE
      )
    } else if (grepl(".tif", target_S)) {
      target_S <- terra::rast(file.path(hydroweight_dir, target_S))
    }
  }

  ## write target_S and adjust to clip_region if sf
  if (class(target_S)[1] == "sf") {
    if (any(sf::st_is(target_S, c("POLYGON", "MULTIPOLYGON")))) {
      target_S_r <- terra::rasterize(sf::as_Spatial(target_S),
                                     y = dem_clip,
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
                                     field = 1,
                                     overwrite = TRUE
      )
    }

    if (any(sf::st_is(target_S, c("LINESTRING", "MULTILINESTRING")))) {
      dem_clip_extent <- sf::st_as_sfc(sf::st_bbox(terra::ext(dem_clip)))
      sf::st_crs(dem_clip_extent) <- sf::st_crs(dem_crs)

      target_S_int <- sf::st_intersects(target_S, dem_clip_extent)
      target_S_int <- target_S[lengths(target_S_int) > 0, ]

      target_S_r <- terra::rasterize(sf::as_Spatial(target_S_int),
                                     y = dem_clip,
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
                                     field = 1,
                                     overwrite = TRUE
      )
    }

    if (sf::st_is(target_S, "POINT")) {
      target_S_r <- terra::rasterize(sf::st_coordinates(target_S),
                                     y = dem_clip,
                                     field = 1,
                                     filename = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
                                     overwrite = TRUE
      )
    }
  }

  ## write target_S and adjust to clip_region if RasterLayer
  if (class(target_S)[1] %in% c("SpatRaster","RasterLayer")) {
    if (class(target_S)[1] %in% c("RasterLayer")) target_S_r <- terra::rast(target_S)
    if (class(target_S)[1] %in% c("SpatRaster")) target_S_r <- target_S
    target_S_r <- terra::project(target_S_r, dem_clip, method = "near",overwrite = TRUE,filename=file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")))
    terra::writeRaster(target_S_r, file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
                       overwrite = TRUE
    )
  }

  ## Prepare OS_combine ----
  if (is.null(OS_combine)) {
    OS_combine <- FALSE
  }

  if (OS_combine == TRUE) {
    (target_S_OS <- target_S_r)
    (target_O_OS <- target_O_r)

    target_S_OS[target_S_OS > 0] <- 1
    target_O_OS[target_O_OS > 0] <- 1

    mosaic_list <- list(target_O_OS, target_S_OS)
    mosaic_list$fun <- sum
    mosaic_list$na.rm <- TRUE

    OS_combine_r <- do.call(terra::mosaic, mosaic_list)
    OS_combine_r[OS_combine_r > 0] <- 1
    OS_combine_r[OS_combine_r == 0] <- NA

    terra::crs(OS_combine_r) <- dem_crs
    terra::writeRaster(OS_combine_r,
                       file.path(hydroweight_dir, paste0(target_uid,"_TEMP-OS_combine.tif")),
                       overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                       NAflag = -9999
    )
  }

  ## RUN DISTANCE-WEIGHTING ----

  message("Running distance-weighting @ ", Sys.time())

  if ("lumped" %in% weighting_scheme) {
    # whitebox::wbt_reclass(
    #   input = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")),
    #   output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-lumped.tif")),
    #   reclass_vals = "1,min,1000000"
    # ) # because max sometimes misses max values?

    # lumped_inv <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-lumped.tif")))
    # terra::crs(lumped_inv) <- dem_crs
    # lumped_inv <- terra::setValues(terra::rast(lumped_inv), lumped_inv[])

    lumped_inv <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")))
    lumped_inv[!is.na(lumped_inv)]<-1
  }

  ## iEucO, Euclidean distance to target_O ----
  if ("iEucO" %in% weighting_scheme) {
    # whitebox::wbt_reclass(
    #   input = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")),
    #   output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
    #   reclass_vals = "1,min,1000000"
    # ) # because max sometimes misses max values?

    lumped_inv <- dem_clip
    lumped_inv[!is.na(lumped_inv)]<-1
    terra::writeRaster(lumped_inv,file.path(hydroweight_dir, paste0(target_uid,"_TEMP_dem_clip_cost.tif"),overwrite=T))

    whitebox::wbt_cost_distance(
      source = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
      cost = file.path(hydroweight_dir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
      out_accum = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")),
      out_backlink = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_backlink.tif")),
      verbose_mode = FALSE
    )

    iEucO <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")))
    terra::crs(iEucO) <- dem_crs

    if(is.list(inv_function)){

      if("iEucO" %in% names(inv_function)){

        iEucO_inv <- terra::app(iEucO, fun = inv_function$iEucO)

      } else {

        stop("no `iEucO` found in names(inv_function); please name elements of `inv_function`")

      }


    } else {

      iEucO_inv <- terra::app(iEucO, fun = inv_function)

    }

    terra::crs(iEucO_inv) <- dem_crs
    terra::writeRaster(iEucO_inv,
                       file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iEucO.tif")),
                       overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                       NAflag = -9999
    )
  }

  ## iEucS, Euclidean distance to target_S ----
  if ("iEucS" %in% weighting_scheme) {
    if (OS_combine == FALSE) {
      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
        cost = file.path(hydroweight_dir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
        out_accum = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")),
        out_backlink = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_backlink.tif")),
        verbose_mode = FALSE
      )

      iEucS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")))
      terra::crs(iEucS) <- dem_crs

      if(is.list(inv_function)){

        if("iEucS" %in% names(inv_function)){

          iEucS_inv <- terra::app(iEucS, fun = inv_function$iEucS)

        } else {

          stop("no `iEucS` found in names(inv_function); please name elements of `inv_function`")

        }


      } else {

        iEucS_inv <- terra::app(iEucS, fun = inv_function)

      }

      terra::crs(iEucS_inv) <- dem_crs
      terra::writeRaster(iEucS_inv,
                         file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iEucS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if (OS_combine == TRUE) {
      whitebox::wbt_cost_distance(
        source = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-OS_combine.tif")),
        cost = file.path(hydroweight_dir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
        out_accum = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")),
        out_backlink = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_backlink.tif")),
        verbose_mode = FALSE
      )

      iEucS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-cost_distance.tif")))
      terra::crs(iEucS) <- dem_crs

      if(is.list(inv_function)){

        if("iEucS" %in% names(inv_function)){

          iEucS_inv <- terra::app(iEucS, fun = inv_function$iEucS)

        } else {

          stop("no `iEucS` found in names(inv_function); please name elements of `inv_function`")

        }


      } else {

        iEucS_inv <- terra::app(iEucS, fun = inv_function)

      }

      terra::crs(iEucS_inv) <- dem_crs
      terra::writeRaster(iEucS_inv,
                         file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iEucS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }
  }

  ## iFLO, flow line distance to target_O ----
  if ("iFLO" %in% weighting_scheme | "HAiFLO" %in% weighting_scheme) {
    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")),
      streams = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_O_clip.tif")),
      output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLO.tif")),
      verbose_mode = FALSE
    )

    iFLO <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLO.tif")))
    terra::crs(iFLO) <- dem_crs

    if(is.list(inv_function)){

      if("iFLO" %in% names(inv_function)){

        iFLO_inv <- terra::app(iFLO, fun = inv_function$iFLO)

      } else {

        stop("no `iFLO` found in names(inv_function); please name elements of `inv_function`")

      }


    } else {

      iFLO_inv <- terra::app(iFLO, fun = inv_function)

    }

    terra::crs(iFLO_inv) <- dem_crs
    terra::writeRaster(iFLO_inv,
                       file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iFLO.tif")),
                       overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                       NAflag = -9999
    )
  }

  ## iFLS, flow line distance to target_S ----
  if ("iFLS" %in% weighting_scheme | "HAiFLS" %in% weighting_scheme) {
    if (OS_combine == FALSE) {
      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")),
        streams = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-target_S_clip.tif")),
        output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")),
        verbose_mode = FALSE
      )

      iFLS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")))
      terra::crs(iFLS) <- dem_crs

      if(is.list(inv_function)){

        if("iFLS" %in% names(inv_function)){

          iFLS_inv <- terra::app(iFLS, fun = inv_function$iFLS)

        } else {

          stop("no `iFLS` found in names(inv_function); please name elements of `inv_function`")

        }


      } else {

        iFLS_inv <- terra::app(iFLS, fun = inv_function)

      }

      terra::crs(iFLS_inv) <- dem_crs
      terra::writeRaster(iFLS_inv,
                         file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iFLS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if (OS_combine == TRUE) {
      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-dem_clip.tif")),
        streams = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-OS_combine.tif")),
        output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")),
        verbose_mode = FALSE
      )

      iFLS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")))
      terra::crs(iFLS) <- dem_crs

      if(is.list(inv_function)){

        if("iFLS" %in% names(inv_function)){

          iFLS_inv <- terra::app(iFLS, fun = inv_function$iFLS)

        } else {

          stop("no `iFLS` found in names(inv_function); please name elements of `inv_function`")

        }


      } else {

        iFLS_inv <- terra::app(iFLS, fun = inv_function)

      }

      terra::crs(iFLS_inv) <- dem_crs
      terra::writeRaster(iFLS_inv,
                         file.path(hydroweight_dir, paste0(target_uid,"_TEMP-iFLS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }
  }

  ## HA-iFLO and HA-iFLS, hydrologically active ----
  if (("HAiFLO" %in% weighting_scheme) | ("HAiFLS" %in% weighting_scheme)) {

    if(is.null(flow_accum)){
      stop("`flow_accum` required for HAiFLO and/or HAiFLS calculations")
    }

    flow_accum_clip_path<-file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flow_accum_clip.tif"))


    flow_accum_clip<-terra::crop(
      x=terra::rast(file.path(hydroweight_dir, flow_accum)),
      y=terra::vect(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.shp")))
    )
    terra::crs(flow_accum_clip)<-dem_crs
    flow_accum_clip<-terra::project(flow_accum_clip,dem_clip, method = "near",overwrite = TRUE,filename=file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flow_accum_clip.tif")))
    terra::writeRaster(flow_accum_clip,
                       flow_accum_clip_path,
                       overwrite=T)


    # whitebox::wbt_clip_raster_to_polygon( #PS I think this is slower than terra
    #   input = file.path(hydroweight_dir, flow_accum),
    #   polygons = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-clip_region.shp")),
    #   output = file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flow_accum_clip.tif")),
    #   verbose_mode = FALSE
    # )

    accum_clip <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flow_accum_clip.tif")))
    # accum_clip <- accum_clip + 1 # PS: I don't know if this is necessary anymore
    terra::crs(accum_clip) <- dem_crs

    accum_clip<-terra::project(accum_clip,dem_clip,method="near",overwrite = TRUE,filename=file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLO.tif")))

    if ("HAiFLO" %in% weighting_scheme) {

      HAiFLO <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLO.tif")))
      terra::crs(HAiFLO) <- dem_crs

      if(is.list(inv_function)){

        if("HAiFLO" %in% names(inv_function)){

          iFLO_inv <- terra::app(HAiFLO, fun = inv_function$HAiFLO)

        } else {

          stop("no `HAiFLO` found in names(inv_function); please name elements of `inv_function`")

        }


      } else {

        HAiFLO_inv <- terra::app(HAiFLO, fun = inv_function)

      }

      HAiFLO_inv <- HAiFLO_inv * accum_clip

      terra::crs(HAiFLO_inv) <- dem_crs
      terra::writeRaster(HAiFLO_inv,
                         file.path(hydroweight_dir, paste0(target_uid,"_TEMP-HAiFLO.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if ("HAiFLS" %in% weighting_scheme) {
      if (OS_combine == TRUE) {

        HAiFLS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")))
        terra::crs(HAiFLS) <- dem_crs

        if(is.list(inv_function)){

          if("HAiFLS" %in% names(inv_function)){

            iFLS_inv <- terra::app(HAiFLS, fun = inv_function$HAiFLS)

          } else {

            stop("no `HAiFLS` found in names(inv_function); please name elements of `inv_function`")

          }


        } else {

          HAiFLS_inv <- terra::app(HAiFLS, fun = inv_function)

        }

        HAiFLS_inv <- HAiFLS_inv * accum_clip
        HAiFLS_inv <- terra::mask(HAiFLS_inv, OS_combine_r, maskvalues = 1)

        terra::writeRaster(HAiFLS_inv,
                           file.path(hydroweight_dir, paste0(target_uid,"_TEMP-HAiFLS.tif")),
                           overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                           NAflag = -9999
        )
      }

      if (OS_combine == FALSE) {
        HAiFLS <- terra::rast(file.path(hydroweight_dir, paste0(target_uid,"_TEMP-flowdist-iFLS.tif")))
        terra::crs(HAiFLS) <- dem_crs

        if(is.list(inv_function)){

          if("HAiFLS" %in% names(inv_function)){

            iFLS_inv <- terra::app(HAiFLS, fun = inv_function$HAiFLS)

          } else {

            stop("no `HAiFLS` found in names(inv_function); please name elements of `inv_function`")

          }


        } else {

          HAiFLS_inv <- terra::app(HAiFLS, fun = inv_function)

        }

        HAiFLS_inv <- HAiFLS_inv * accum_clip
        HAiFLS_inv <- terra::mask(HAiFLS_inv, target_S_r, maskvalues = 1)

        terra::writeRaster(HAiFLS_inv,
                           file.path(hydroweight_dir, paste0(target_uid,"_TEMP-HAiFLS.tif")),
                           overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                           NAflag = -9999
        )
      }
    }
  }

  ## PREPARE OUTPUT ----

  weighting_scheme_inv <- paste0(weighting_scheme, "_inv")

  dist_list <- lapply(weighting_scheme_inv, function(x) {
    get(x) #I think this can cause problems
  })
  names(dist_list) <- weighting_scheme

  dist_list<-lapply(dist_list,terra::wrap) # will need to use wrap() here for terra: https://github.com/rspatial/terra/issues/50

  saveRDS(dist_list, file = file.path(
    hydroweight_dir,
    paste0(target_uid, "_inv_distances.rds")
  ))

  ## CLEAN UP ----
  temp_rasters <- list.files(path = hydroweight_dir, pattern = paste0(target_uid,"_hydroweight_"),
                             full.names = TRUE)
  file.remove(temp_rasters)

  temp_rasters <- list.files(path = hydroweight_dir, pattern = paste0(target_uid,"_TEMP-"),
                             full.names = TRUE)
  file.remove(temp_rasters)

  return(dist_list)
}
