#' Align spatial data for `hydroweight::` functions
#'
#' `hydroweight::process_input()` is an internal function that aligns spatial data into `SpatVector` or `SpatRaster` types depending on `align_to` class.
#'
#' `input` data is coerced to type `SpatVector` or `SpatRaster` depending on class of `align_to`:
#'    * If `input` is a raster, `align_to` is a raster, and `resample_type` == `bilinear` (e.g., for continuous data), then class `SpatRaster` is output with rasterizing and reprojection occurring as necessary;
#'    * If `input` is a raster, `align_to` is a raster, and `resample_type` == `near` (e.g., for categorical data), then unique values of `input` are coerced into separate binary layers of a `SpatRaster`, with 1 indicating the category/value is present;
#'    * If `input` is a vector, `align_to` is a raster, unique values of `input_variable_names` are each coerced into separate layers in a `SpatRaster`, with 1 indicating the category/value is present.

#'    , and `align_to` is of class `SpatRaster`, `input` is coerced to class `SpatRaster`, rasterizing and/or reprojecting as necessary;
#'    * If `resample_type` == `near`, (e.g., for categorical data), `align_to` is of class `SpatRaster`, and `input` is of class `SpatRaster`,  and
#'    * If `input` is of class `SpatVector` and `align_to` is of class `SpatRaster`,
#'    * If `input` is of class `SpatVector` and `align_to` is
#'
#'  `snap` can have very important implications for hydroweight and
#'
#' @param input Input spatial data: `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/input.shp").
#' @param input_name Input name only used in error or warning messages: `NULL` or `character`.
#' @param input_variable_names `input` variable names that should be included in output: `NULL` or `character`. If `NULL`, all variables in `input` are used.
#' @param align_to Spatial data that `input` will be aligned to: `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster` or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/align_to.tif"). If `NULL`, all variables in `input` are returned as `SpatVector` or `SpatRaster`.
#' @param clip_region `input` will be clipped and masked to this region: `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or `character` (full file path with extension, e.g., "C:/Users/Administrator/Desktop/clip_region.shp"). Internally converted to `SpatVector`. If `NULL`, full extent of `input` is returned.
#' @param resample_type Type of resampling needed to conver to raster (typically categorical = `near` and numerical = `bilinear`): `character`. If `bilinear` (default), the `input_variable_names` being coerced are numeric, if `'near'`, the `input_variable_names` being coerced are categorical, for the purposes of resampling/projecting if `align_to` is `SpatRaster`.
#' @param snap Type of snap scheme to use during `terra::crop()`: `character`. One of "near", "in", or "out". Used to align y to the geometry of x. Default is `near`.
#' @param working_dir Folder path for temporary file storage: `NULL` or `character`. If `NULL`, assigned internally.
#' @param ... other variables passed to various `terra` functions.
#' @return an object of class `SpatVector` or `SpatRaster` depending on `align_to`
#' @export

process_input <- function(input = NULL,
                          input_name = NULL,
                          input_variable_names = NULL,
                          align_to = NULL,
                          clip_region = NULL,
                          resample_type = c("bilinear", "near"),
                          snap = c("near", "in", "out"),
                          working_dir = NULL,
                          ...) {

  ## Set up outputs and working directories ------------------------------------

  ## Output
  if (is.null(input)) {
    return(input)
  }
  output <- NULL

  ## Working directories
  if (is.null(working_dir)) {
    tdir <- file.path(gsub("file", "", tempfile()))
    working_dir <- tdir # also make working_dir
  } else {
    tdir <- file.path(working_dir, basename(gsub("file", "", tempfile())))
  }
  if (!dir.exists(tdir)) {
    dir.create(tdir)
  }

  ## Set terra options
  terra::terraOptions(tempdir = tdir, verbose = FALSE)

  ## match.arg
  resample_type <- match.arg(resample_type)
  snap <- match.arg(snap)

  # Reconcile Types ---------------------------------------------------------
  if (inherits(input, "character")) {
    if (grepl("\\.shp$", input)) {
      output <- terra::vect(input) # terra seems to have trouble reading crs sometimes (sf::read_sf())
    } else if (grepl("\\.tif$|\\.tiff", input)) {
      output <- terra::rast(input)
    } else {
      stop(paste0(input_name, " must be spcified as vector or raster layer, or a character string ending in .shp or .tif"))
    }
  }
  if (inherits(input, "RasterLayer")) {
    output <- terra::rast(input)
  }
  if (inherits(input, "PackedSpatRaster")) {
    output <- terra::rast(input)
  }
  if (inherits(input, "PackedSpatVector")) {
    output <- terra::vect(input)
  }
  if (inherits(input, c("sf", "sfc"))) {
    output <- terra::vect(input)
  }
  if (inherits(input, c("SpatRaster", "SpatVector"))) {
    output <- input
  }
  if (is.na(terra::crs(output)) | is.null(terra::crs(output))) {
    stop("'output' crs() is NULL or NA. Apply a CRS before continuing")
  }

  ## Construct output using input_variable_names as necessary
  orig_output <- output

  if (is.null(input_variable_names)) {
    input_variable_names <- names(output)
  }
  if (any(!input_variable_names %in% names(output))) {
    stop("some 'input_variable_names' not in input")
  }
  if (length(input_variable_names) > 0) {
    if (inherits(output, "SpatVector")) {
      output <- output[, input_variable_names]
    }
    if (inherits(output, "SpatRaster")) {
      output <- terra::subset(output, input_variable_names)
    }
  }

  ## Align to align_to ---------------------------------------------------------

  if (!is.null(align_to)) {
    ## This seems quite circular, given that it is within the process_input
    ## function, check in on whether something different needs to be done.

    align_to <- process_input(input = align_to, working_dir = tdir) # builds a temp directory

    clip_region <- process_input(clip_region,
      align_to = terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",
      crs = terra::crs(align_to)
      ),
      working_dir = tdir
    )

    if (is.na(terra::crs(align_to)) | is.null(terra::crs(align_to))) {
      stop("'align_to' crs() is NULL or NA. Apply CRS before continuing")
    }

    ## If align_to is raster ---------------------------------------------------
    if (inherits(align_to, "SpatRaster")) {

      ## If output (i.e., same class as input) is
      if (inherits(output, "SpatVector")) {

        output <- terra::project(output, align_to)

        ## For categorical vector inputs ---------------------------------------
        if (resample_type == "near") {
          # output_split <- lapply(stats::setNames(input_variable_names, input_variable_names), function(x) {
          #   out <- output %>%
          #     sf::st_as_sf() %>%
          #     dplyr::select(tidyselect::any_of(x)) %>%
          #     split(.[[x]]) %>%
          #     lapply(terra::vect)
          #
          #   fl <- lapply(out, function(x) {
          #     file.path(tdir, paste0(basename(tempfile()), ".shp"))
          #     })
          #   sv <- lapply(names(out), function(x) {
          #     terra::writeVector(out[[x]], filename = fl[[x]], overwrite = T)
          #     })
          #
          #   out <- lapply(fl, terra::vect)
          #
          #   out <- lapply(out, function(y) {
          #     names(y) <- paste0(x, "_", unlist(y[[1]])[[1]])
          #     return(y)
          #   })
          #
          #   return(out)
          # })

          ## More package safe version of code above
          output_split <- lapply(stats::setNames(
            input_variable_names,
            input_variable_names
          ), function(x) {
            # Convert to sf
            sfobj <- sf::st_as_sf(output)

            # Select only the column(s) of interest
            sfobj <- sfobj[, x, drop = FALSE]

            # Split by the grouping column (x)
            out <- split(sfobj, sfobj[[x]]) # <- no dot pronoun here
            out <- lapply(out, terra::vect)

            # Create file paths
            fl <- lapply(out, function(y) {
              file.path(tdir, paste0(basename(tempfile()), ".shp"))
            })

            # Write each group to file
            sv <- lapply(names(out), function(nm) {
              terra::writeVector(out[[nm]], filename = fl[[nm]], overwrite = TRUE)
            })

            # Reload as SpatVector
            out <- lapply(fl, terra::vect)

            # Rename layers
            out <- lapply(out, function(y) {
              names(y) <- paste0(x, "_", unlist(y[[1]])[[1]])
              return(y)
            })

            return(out)
          })

          output_split_nms <- sapply(names(output_split), function(x) {
            paste0(x, "_", names(output_split[[x]]))
          })

          output_split <- unlist(output_split, recursive = F)
          names(output_split) <- output_split_nms

          output <- lapply(output_split, function(x) {
            fl <- file.path(tdir, paste0(basename(tempfile()), ".tif"))
            out <- terra::rasterize(
              x = x,
              y = align_to,
              field = "",
              overwrite = TRUE,
              ...
            )

            sv <- terra::writeRaster(out, filename = fl, overwrite = TRUE, gdal = "COMPRESS=NONE")

            return(terra::rast(fl))
          })

          output <- lapply(stats::setNames(names(output), names(output)), function(x) {
            out <- output[[x]]
            names(out) <- x
            return(out)
          })
        }

        ## For numeric vector inputs -------------------------------------------
        if (resample_type == "bilinear") {
          output <- lapply(input_variable_names, function(x) {
            fl <- file.path(tdir, paste0(basename(tempfile()), ".tif"))
            out <- terra::rasterize(
              x = output, # PS: Development version of terra gives a warning here, may cause problems in the future
              y = align_to,
              field = x,
              overwrite = TRUE,
              ...
            )
            sv <- terra::writeRaster(out, filename = fl, overwrite = TRUE, gdal = "COMPRESS=NONE")
            return(terra::rast(fl))
          })
        }

        output <- terra::rast(output)
      }
    }

    ## If align_to is vector ---------------------------------------------------
    if (inherits(align_to, "SpatVector")) {
      if (inherits(output, "SpatRaster")) {
        if (terra::is.polygons(align_to)) {
          output <- terra::as.polygons(output, ...)
        }
        if (terra::is.points(align_to)) {
          output <- terra::as.points(output, ...)
        }
        if (terra::is.lines(align_to)) {
          output <- terra::as.lines(output, ...)
        }
      }
    }

    ## Reprojecting to align_to as necessary -----------------------------------

    ## Compare output and aligh_to
    if (inherits(output, "SpatRaster")) {
      need_reproj <- terra::compareGeom(
        x = output,
        y = align_to,
        lyrs = FALSE,
        crs = TRUE,
        ext = TRUE,
        rowcol = TRUE,
        res = TRUE,
        warncrs = FALSE,
        stopOnError = FALSE,
        messages = FALSE
      )
      need_reproj <- !need_reproj # take inverse

      ## Reproject as necessary
      if (need_reproj) {
        if (length(resample_type) > 1) { # BK: maybe this should be before doing anything?
          stop("'resample_type' must be one of: 'bilinear','near'")
        }

        output <- terra::project(
          x = output,
          y = align_to,
          method = resample_type,
          overwrite = TRUE,
          ...
        )
      }
    }

    if (inherits(output, "SpatVector")) {
      output <- terra::project(
        x = output,
        y = align_to,
        ...
      )
    }
  }

  ## Clip and mask to clip_region -----------------------------------------------

  if (!is.null(clip_region) & !is.null(align_to)) {
    ## For vector data
    if (inherits(output, "SpatVector")) {
      output <- terra::crop(
        x = output,
        y = clip_region
      )
    }

    ## For raster data
    if (inherits(output, "SpatRaster")) {
      output <- terra::crop(
        x = output,
        snap = snap,
        y = clip_region,
        mask = TRUE,
        overwrite = TRUE
      )

      output <- terra::mask(
        x = output,
        mask = clip_region,
        overwrite = TRUE
      )
    }
  }

  ## Split output into separate layers if multiple numerical values present ----
  if (inherits(output, "SpatRaster")) {
    if (resample_type == "near") {
      if (terra::nlyr(output) > 1) {
        brick_list <- terra::split(output, names(output))
        names(brick_list) <- names(output)
      } else {
        brick_list <- output
      }

      brick_list <- lapply(names(brick_list), function(y) {
        uv <- unlist(terra::unique(brick_list[[y]]))
        uv <- uv[!sapply(uv, is.nan)]
        uv <- uv[!sapply(uv, is.na)]
        uv <- uv[!sapply(uv, is.infinite)]

        if (length(uv) > 1) {
          out <- lapply(uv, function(x) {
            repl <- terra::values(brick_list[[y]])
            repl[repl != x] <- NA
            repl[!is.na(repl)] <- 1
            out <- terra::setValues(brick_list[[y]], repl)
            if (length(uv) > 1) {
              names(out) <- paste0(y, "_", x)
            }
            if (length(uv) == 1) {
              names(out) <- y
            }
            return(out)
          })

          nms <- sapply(out, names)
          out <- terra::rast(out)
          names(out) <- nms
        } else {
          out <- brick_list[[y]]
        }
        return(out)
      })

      output <- terra::rast(brick_list)
    }
  }

  if (inherits(output, "SpatRaster")) {
    if (!inherits(orig_output, "SpatVector")) {
      need_save <- terra::compareGeom(
        x = output,
        y = orig_output,
        lyrs = T,
        crs = T,
        ext = T,
        rowcol = T,
        res = T,
        warncrs = F,
        stopOnError = F,
        messages = F
      )
      need_save <- !need_save # take inverse
    } else {
      need_save <- TRUE
    }

    if (need_save) {
      final_temp <- file.path(tdir, paste0(basename(tempfile()), ".tif"))
      terra::writeRaster(output, final_temp, overwrite = TRUE, gdal = "COMPRESS=NONE")
    }
  }

  if (inherits(output, "SpatVector")) {
    final_temp <- file.path(tdir, paste0(basename(tempfile()), ".shp"))
    terra::writeVector(output, final_temp, overwrite = T)
  }

  return(output)
}
