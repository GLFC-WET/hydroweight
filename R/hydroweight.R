#' Generate inverse distance-weighted rasters
#'
#' `hydroweight::hydroweight()` generates distance-weighted rasters for targets
#' on a digital elevation model raster. Examples of targets include single points,
#' areas such as lakes, or linear features such as streams. The function outputs
#' a list of `weighting_scheme` and/or an accompanying `.zip` file
#' of distance-weighted rasters for targets (`target_O` is a point/area target as
#' in iFLO and `target_S` is a stream/linear feature target as in iFLS in
#' Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function acts on a single set of targets but can produce
#' multiple weights. The distance-weighted rasters, can be used for generating
#' distance-weighted landscape statistics using `hydroweight::hydroweight_attributes()` (e.g., % urban
#' cover weighted by flow distance to a point). See https://github.com/GLFC-WET/hydroweight for workflows.
#'
#' Spatial layers should align (i.e., identical coordinate reference systems (CRS) and identical dem and flow_accum extents, resolutions, and CRS).
#' Through processing, targets are converted to the resolution and CRS of the
#' digital elevation model/flow accumulation CRS. _TEMP_* files are generated in
#'  `hydroweight_dir` depending on processing step and can be overwritten.
#'
#' For `weighting_scheme`:
#'
#'  * "lumped" indicates all weights = 1
#'
#'  * "iEucO" indicates Euclidean distance to `target_O`
#'
#'  * "iEucS" indicates Euclidean distance to `target_S`
#'
#'  * "iFLO" indicates d8 flow-path distance to `target_O`
#'
#'  * "iFLS" indicates d8 flow-path distance to `target_S`
#'
#'  * "HAiFLO" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to `target_O`
#'
#'  * "HAiFLS" indicates d8 hydrologically-active (proportional to flow accumulation) flow-path distance to `target_S`
#'
#' @param hydroweight_dir File path to write resulting `.zip` file: `character`.
#' @param target_O Target for iEucO or iFLO: `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster` or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/sample_point.shp").
#' @param target_S Target for iEucS or iFLS: `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster` or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/sample_streams.shp").
#' @param target_uid Unique identifier to precede intermediate/exported (i.e., "target_uid"_inv_distances.zip): `character`.
#' @param OS_combine Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use `TRUE` or `FALSE`. If `TRUE`, this allows cells surrounding `target_O` to flow directly into `target_O` rather than be forced through `target_S`.
#' @param clip_region Spatial inputs will be clipped to this region: `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/clip_region.shp"). Internally converted to `SpatVector`. If `NULL`, full extent of inputs are returned.
#' @param dem Digital elevation model raster: `RasterLayer`, `SpatRaster`, `PackedSpatRaster` or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif").
#' @param flow_accum  Flow accumulation raster (units: # of cells): `RasterLayer`, `SpatRaster`, `PackedSpatRaster` or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/flow_accum.tif"). Should align with `dem`.
#' @param weighting_scheme One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS").
#' @param inv_function Inverse function used in `terra:app()` to convert distances to inverse distances. Format is function or named list of functions based on `weighting_scheme` names.  Default: `(X * 0.001 + 1)^-1` assumes projection is in distance units of m and converts to distance units of km.
#' @param snap Type of snap scheme to use during `terra::crop()`: `character`. One of "near", "in", or "out". Used to align y to the geometry of x. Default is `near`.
#' @param return_products If `TRUE`, a list containing the file path to write resulting `.zip` file, and `distance_weights` rasters. If `FALSE`, file path only.
#' @param save_output Should output rasters be saved to a zip file?
#' @param wrap_return_products Should `return_products` be `terra::wrap()` (necessary for running function in parallel)? If `TRUE`, wrapped.
#' @param clean_tempfiles Should temporary files be removed? If `TRUE`, temporary files are removed.
#'
#' @return Named list of `SpatRaster` or `PackedSpatRaster` (if `wrap_return_products` = `TRUE`) distance-weighted rasters and location of accompanying \code{*.zip} in \code{hydroweight_dir}
#' @export
#'
hydroweight <- function(
    hydroweight_dir = NULL,
    target_O = NULL,
    target_S = NULL,
    target_uid = NULL,
    OS_combine = NULL,
    clip_region = NULL,
    dem,
    flow_accum = NULL,
    weighting_scheme = c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
    inv_function = function(x) { (x * 0.001 + 1)^-1 },
    snap = c("near", "in", "out"),
    return_products = TRUE,
    wrap_return_products = TRUE,
    save_output = TRUE,
    clean_tempfiles = TRUE
) {

  ## SETUP ---------------------------------------------------------------------

  message("Preparing hydroweight layers @ ", Sys.time())

  ## Initial checks
  if (is.null(target_uid)) {
    stop("`target_uid` must be specified")
  }

  if (!save_output && !return_products) {
    stop("'return_products' and 'save_output' cannot both be FALSE")
  }

  ## Match arguments
  weighting_scheme <- match.arg(weighting_scheme, several.ok = TRUE)
  snap <- match.arg(snap)

  ## Use a stable unique working subdir
  if (!is.null(hydroweight_dir)) {
    own_tempdir <- file.path(hydroweight_dir, paste0(target_uid, "_hw_tmp_", basename(tempfile())))
  } else {
    own_tempdir <- file.path(tempdir(), paste0(target_uid, "_hw_tmp_", basename(tempfile())))
  }
  dir.create(own_tempdir, recursive = TRUE, showWarnings = FALSE)

  ## Set terra tempdir for out-of-memory rasters
  old_terra_opts <- terra::terraOptions()
  terra::terraOptions(tempdir = own_tempdir, verbose = FALSE)

  ## Ensure cleanup / terra options restoration happen even on error
  on.exit({
    # restore terra options
    do.call(terra::terraOptions, old_terra_opts)
    if (isTRUE(clean_tempfiles)) {
      unlink(own_tempdir, recursive = TRUE, force = TRUE)
      suppressWarnings(terra::tmpFiles(current = TRUE, orphan = TRUE, old = TRUE, remove = TRUE))
    }
  }, add = TRUE)

  ## Whitebox options (quiet)
  whitebox::wbt_options(verbose = FALSE)

  ## Generate inverse function helper to save downwstream code
  apply_inverse <- function(r, key) {
    if (is.list(inv_function)) {
      if (!key %in% names(inv_function)) {
        stop("No `", key, "` found in names(inv_function); please name elements of `inv_function` accordingly.")
      }
      terra::app(r, fun = inv_function[[key]])
    } else {
      terra::app(r, fun = inv_function)
    }
  }

  ## Generate write raster helper to save downstream code
  write_raster_if <- function(r, filename) {
    terra::writeRaster(r, file.path(own_tempdir, filename),
                       overwrite = TRUE, gdal = "COMPRESS=NONE", NAflag = -9999)
  }

  ## PROCESS INPUTS ------------------------------------------------------------

  ## dem
  dem <- process_input(input = dem, align_to = dem, input_name = "dem", working_dir = own_tempdir)
  dem <- dem[[1]]
  dem_crs <- terra::crs(dem)

  ## clip_region
  clip_region <- process_input(
    input        = clip_region,
    align_to     = terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))", crs = dem_crs),
    clip_region  = terra::as.polygons(terra::ext(dem), crs = terra::crs(dem)),
    snap         = snap,
    input_name   = "clip_region",
    working_dir  = own_tempdir
  )
  clip_region <- clip_region[, 1]

  ## flow_accum
  flow_accum <- process_input(
    input         = flow_accum,
    align_to      = dem,
    clip_region   = clip_region,
    resample_type = "bilinear",
    input_name    = "flow_accum",
    working_dir   = own_tempdir
  )
  flow_accum <- flow_accum[[1]]

  ## target_O
  target_O <- process_input(
    input        = target_O,
    align_to     = dem,
    clip_region  = clip_region,
    input_name   = "target_O",
    working_dir  = own_tempdir
  )
  target_O <- target_O[[1]]
  target_O[!is.na(target_O)] <- 1

  ## target_S
  target_S <- process_input(
    input        = target_S,
    align_to     = dem,
    clip_region  = clip_region,
    input_name   = "target_S",
    working_dir  = own_tempdir
  )
  target_S <- target_S[[1]]
  target_S[!is.na(target_S)] <- 1

  ## OS combine raster
  if (OS_combine && !is.null(target_O) && !is.null(target_S)) {
    target_S_OS <- target_S
    target_O_OS <- target_O
    target_S_OS[target_S_OS > 0] <- 1
    target_O_OS[target_O_OS > 0] <- 1
    OS_combine_r <- terra::mosaic(x = terra::sprc(list(target_O_OS, target_S_OS)), fun = "max")
    OS_combine_r[OS_combine_r == 0] <- NA
  } else {
    OS_combine_r <- NULL
  }

  ## WRITE TEMP RASTERS --------------------------------------------------------

  rast_list <- list(
    dem         = dem,
    flow_accum  = flow_accum,
    clip_region = clip_region,
    target_O    = target_O,
    target_S    = target_S,
    OS_combine  = OS_combine_r
  )

  rast_names <- c(
    paste0(target_uid, "_TEMP_dem_clip.tif"),
    paste0(target_uid, "_TEMP_flow_accum_clip.tif"),
    paste0(target_uid, "_TEMP_clip_region.tif"),
    paste0(target_uid, "_TEMP_target_O_clip.tif"),
    paste0(target_uid, "_TEMP_target_S_clip.tif"),
    paste0(target_uid, "_TEMP_OS_combine.tif")
  )
  names(rast_list) <- rast_names

  # If a list element is a vector, write as .shp, else .tif
  conv_names <- lapply(names(rast_list), function(x) {
    if (inherits(rast_list[[x]], "SpatVector")) {
      gsub("\\.tif$", ".shp", x)
    } else {
      x
    }
  })
  names(rast_list) <- conv_names

  # Drop NULLs and zero-length
  rast_list <- rast_list[!sapply(rast_list, is.null)]
  rast_list <- rast_list[sapply(rast_list, length) > 0]

  for (nm in names(rast_list)) {
    if (inherits(rast_list[[nm]], "SpatVector")) {
      terra::writeVector(rast_list[[nm]], file.path(own_tempdir, nm), overwrite = TRUE)
    } else {
      terra::writeRaster(rast_list[[nm]], file.path(own_tempdir, nm),
                         overwrite = TRUE, gdal = "COMPRESS=NONE")
    }
  }

  ## RUN DISTANCE-WEIGHTING ----------------------------------------------------
  message("Running distance-weighting @ ", Sys.time())

  ## Generate lumped -----------------------------------------------------------
  lumped_inv <- dem
  lumped_inv[!is.na(lumped_inv)] <- 1
  names(lumped_inv) <- "lumped"
  terra::writeRaster(lumped_inv,
                     file.path(own_tempdir, paste0(target_uid, "_TEMP_dem_clip_cost.tif")),
                     overwrite = TRUE, gdal = "COMPRESS=NONE")
  if (save_output) {
    terra::writeRaster(lumped_inv,
                       file.path(own_tempdir, paste0(target_uid, "_TEMP_lumped.tif")),
                       overwrite = TRUE, gdal = "COMPRESS=NONE")
  }

  ## Generate iEucO ------------------------------------------------------------
  if ("iEucO" %in% weighting_scheme) {
    whitebox::wbt_cost_distance(
      source       = file.path(own_tempdir, paste0(target_uid, "_TEMP_target_O_clip.tif")),
      cost         = file.path(own_tempdir, paste0(target_uid, "_TEMP_dem_clip_cost.tif")),
      out_accum    = file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_distance.tif")),
      out_backlink = file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_backlink.tif")),
      verbose_mode = FALSE
    )
    iEucO <- terra::rast(file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_distance.tif")))
    terra::crs(iEucO) <- dem_crs

    iEucO_inv <- apply_inverse(iEucO, "iEucO")
    terra::crs(iEucO_inv) <- dem_crs
    iEucO_inv <- terra::mask(iEucO_inv, lumped_inv)
    names(iEucO_inv) <- "iEucO"
    if (save_output) write_raster_if(iEucO_inv, paste0(target_uid, "_TEMP_iEucO.tif"))
  }

  ## Generate iEucS ------------------------------------------------------------
  if ("iEucS" %in% weighting_scheme) {
    source_eucS <- if (OS_combine) {
      file.path(own_tempdir, paste0(target_uid, "_TEMP_OS_combine.tif"))
    } else {
      file.path(own_tempdir, paste0(target_uid, "_TEMP_target_S_clip.tif"))
    }

    whitebox::wbt_cost_distance(
      source       = source_eucS,
      cost         = file.path(own_tempdir, paste0(target_uid, "_TEMP_dem_clip_cost.tif")),
      out_accum    = file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_distance.tif")),
      out_backlink = file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_backlink.tif")),
      verbose_mode = FALSE
    )
    iEucS <- terra::rast(file.path(own_tempdir, paste0(target_uid, "_TEMP_cost_distance.tif")))
    terra::crs(iEucS) <- dem_crs

    iEucS_inv <- apply_inverse(iEucS, "iEucS")
    terra::crs(iEucS_inv) <- dem_crs
    iEucS_inv <- terra::mask(iEucS_inv, lumped_inv)
    names(iEucS_inv) <- "iEucS"
    if (save_output) write_raster_if(iEucS_inv, paste0(target_uid, "_TEMP_iEucS.tif"))
  }

  ## Generate iFLO / HAiFLO ----------------------------------------------------
  if (any(c("iFLO", "HAiFLO") %in% weighting_scheme)) {
    whitebox::wbt_downslope_distance_to_stream(
      dem     = file.path(own_tempdir, paste0(target_uid, "_TEMP_dem_clip.tif")),
      streams = file.path(own_tempdir, paste0(target_uid, "_TEMP_target_O_clip.tif")),
      output  = file.path(own_tempdir, paste0(target_uid, "_TEMP_flowdist-iFLO.tif")),
      verbose_mode = FALSE
    )
    iFLO <- terra::rast(file.path(own_tempdir, paste0(target_uid, "_TEMP_flowdist-iFLO.tif")))
    terra::crs(iFLO) <- dem_crs

    if ("iFLO" %in% weighting_scheme) {
      iFLO_inv <- apply_inverse(iFLO, "iFLO")
      terra::crs(iFLO_inv) <- dem_crs
      names(iFLO_inv) <- "iFLO"
      if (save_output) write_raster_if(iFLO_inv, paste0(target_uid, "_TEMP_iFLO.tif"))
    }

    if ("HAiFLO" %in% weighting_scheme) {
      HAiFLO_inv <- apply_inverse(iFLO, "HAiFLO")
      HAiFLO_inv <- HAiFLO_inv * flow_accum
      terra::crs(HAiFLO_inv) <- dem_crs
      names(HAiFLO_inv) <- "HAiFLO"
      if (save_output) write_raster_if(HAiFLO_inv, paste0(target_uid, "_TEMP_HAiFLO.tif"))
    }
  }

  ## Generate iFLS / HAiFLS ----------------------------------------------------
  if (any(c("iFLS", "HAiFLS") %in% weighting_scheme)) {
    streams_src <- if (OS_combine) {
      file.path(own_tempdir, paste0(target_uid, "_TEMP_OS_combine.tif"))
    } else {
      file.path(own_tempdir, paste0(target_uid, "_TEMP_target_S_clip.tif"))
    }

    whitebox::wbt_downslope_distance_to_stream(
      dem     = file.path(own_tempdir, paste0(target_uid, "_TEMP_dem_clip.tif")),
      streams = streams_src,
      output  = file.path(own_tempdir, paste0(target_uid, "_TEMP_flowdist-iFLS.tif")),
      verbose_mode = FALSE
    )

    iFLS <- terra::rast(file.path(own_tempdir, paste0(target_uid, "_TEMP_flowdist-iFLS.tif")))
    terra::crs(iFLS) <- dem_crs

    if ("iFLS" %in% weighting_scheme) {
      iFLS_inv <- apply_inverse(iFLS, "iFLS")
      terra::crs(iFLS_inv) <- dem_crs
      names(iFLS_inv) <- "iFLS"
      if (save_output) write_raster_if(iFLS_inv, paste0(target_uid, "_TEMP_iFLS.tif"))
    }

    if ("HAiFLS" %in% weighting_scheme) {
      HAiFLS_inv <- apply_inverse(iFLS, "HAiFLS")
      HAiFLS_inv <- HAiFLS_inv * flow_accum
      if (OS_combine) {
        HAiFLS_inv <- terra::mask(HAiFLS_inv, OS_combine_r, maskvalues = 1)
      } else {
        HAiFLS_inv <- terra::mask(HAiFLS_inv, target_S, maskvalues = 1)
      }
      names(HAiFLS_inv) <- "HAiFLS"
      if (save_output) write_raster_if(HAiFLS_inv, paste0(target_uid, "_TEMP_HAiFLS.tif"))
    }
  }

  ## PREPARE OUTPUTS -----------------------------------------------------------
  weighting_scheme_inv <- paste0(weighting_scheme, "_inv")
  dist_list <- lapply(weighting_scheme_inv, function(x) {
    out <- get(x)
    names(out) <- gsub("_inv", "", x)
    out
  })
  names(dist_list) <- weighting_scheme

  if (save_output) {
    dist_list_out <- lapply(dist_list, function(x) {
      file.rename(
        from = file.path(own_tempdir, paste0(target_uid, "_TEMP_", names(x), ".tif")),
        to   = file.path(own_tempdir, paste0(names(x), ".tif"))
      )
      file.path(own_tempdir, paste0(names(x), ".tif"))
    })

    out_file <- file.path(hydroweight_dir, paste0(target_uid, "_inv_distances.zip"))
    if (file.exists(out_file)) {
      out_file <- file.path(
        hydroweight_dir,
        paste0(target_uid, "_", basename(tempfile()), "_inv_distances.zip")
      )
      warning(paste0("Target .zip file already exists. Saving as: ", out_file))
    }
    utils::zip(out_file, unlist(dist_list_out), flags = "-r9Xjq")
  } else {
    out_file <- NULL
    dist_list_out <- NULL
  }

  if (return_products) {
    if (wrap_return_products) {
      dist_list <- lapply(dist_list, terra::wrap)
    } else {
      dist_list <- dist_list
    }
  } else {
    dist_list <- out_file
  }

  return(dist_list)
}
