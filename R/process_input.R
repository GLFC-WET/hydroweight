#' Aligns spatial data for `hydroweight` functions
#'
#' `hydroweight::process_input()` internally aligns spatial data into `SpatVector` or
#' `SpatRaster` objects depending on the class of `align_to`. The function
#' normalizes inputs (vectors/rasters or file paths), optionally subsets
#' variables, aligns to a target object, and optionally clips/masks to a
#' region.
#'
#' @details
#' **Coercion & alignment behavior (high-level):**
#'
#' \itemize{
#'   \item If `input` is a raster, `align_to` is a raster, and `resample_type == "bilinear"`
#'   (e.g., continuous data), the result is a `SpatRaster` with rasterize/warp
#'   applied as needed.
#'
#'   \item If `input` is a raster, `align_to` is a raster, and `resample_type == "near"`
#'   (e.g., categorical data), unique values of `input` are coerced into
#'   separate binary `SpatRaster` layers (1 = category present).
#'
#'   \item If `input` is a vector and `align_to` is a raster, unique values of
#'   each variable in `input_variable_names` are coerced into separate binary
#'   `SpatRaster` layers (1 = category present).
#' }
#'
#' **Notes:**
#' \itemize{
#'   \item `resample_type` affects reprojection/resampling when `align_to`
#'   is a `SpatRaster`. Use `"bilinear"` for numeric/continuous data and
#'   `"near"` for categorical data.
#'   \item `snap` controls how `terra::crop()` aligns the raster grid to the
#'   clipping geometry (`"near"`, `"in"`, or `"out"`), which can materially
#'   affect hydrologic analyses.
#' }
#'
#' @param input Input spatial data: one of `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or a file path (`character`) with a recognized extension (e.g., `"C:/path/input.shp"`, `"C:/path/input.tif"`).
#' @param input_name Optional label used in error/warning messages. `NULL` or `character`.
#' @param input_variable_names Names of variables in `input` to include in the output: one of `NULL` or `character`. If `NULL`, all variables are used.
#' @param align_to Spatial target to which `input` will be aligned: one of `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or a file path (`character`) with a recognized extension (e.g., `"C:/path/align_to.shp"`, `"C:/path/align_to.tif"`). If `NULL`, the function returns a `SpatVector` or `SpatRaster` in the original input CRS/geometry (after optional subsetting).
#' @param clip_region Optional region for clipping/masking: one of `NULL`, `sf`, `SpatVector`, `PackedSpatVector`, `RasterLayer`, `SpatRaster`, `PackedSpatRaster`, or a file path (`character`) with a recognized extension (e.g., `"C:/path/clip_region.shp"`, `"C:/path/clip_region.tif"`). Internally converted to `SpatVector`. If `NULL`, the full extent of `input` (after alignment) is returned.
#' @param resample_type Type of resampling when converting/alignment to raster (typically categorical = `"near"`, numeric = `"bilinear"`): `character`, one of`c("bilinear","near")`. If `"bilinear"` (default), `input_variable_names` are treated as numeric; if `"near"`, they are treated as categorical for resampling/projecting when `align_to` is a `SpatRaster`.
#' @param snap Snap scheme for `terra::crop()`: `character`, one of `"near"`, `"in"`, or `"out"`. Controls how the raster grid aligns to the clipping geometry. Default is `"near"`.
#' @param working_dir Folder path for temporary file storage: `NULL` or `character`. If `NULL`, a temporary directory is created internally and may be removed at the end.
#' @param persist_final Should the final `SpatRaster` or `SpatVector` be written to disk inside `working_dir` and returned as a file-backed object? Logical, default `FALSE`. Set to `TRUE` when reproducible on-disk outputs or file-backed terra objects are required.
#' @param ... Additional arguments passed to `terra` operations.
#'
#' @return An object of class `SpatVector` or `SpatRaster`, depending on `align_to`.
#' @export
#'
process_input <- function(
    input = NULL,
    input_name = NULL,
    input_variable_names = NULL,
    align_to = NULL,
    clip_region = NULL,
    resample_type = c("bilinear", "near"),
    snap = c("near", "in", "out"),
    working_dir = NULL,
    persist_final = FALSE,
    ...
) {

  ## SETUP ---------------------------------------------------------------------

  if (is.null(input)) {
    return(input)
  }

  ## Match arguments
  resample_type <- match.arg(resample_type)
  snap <- match.arg(snap)

  ## Use a single stable working directory per call
  working_root <- if (is.null(working_dir)) {
    tempfile("hwt_")
  } else {
    working_dir
  }
  dir.create(working_root, showWarnings = FALSE, recursive = TRUE)

  ## Helper for generating temporary file paths
  tmp_file <- function(ext = ".tif") {
    file.path(working_root, paste0("tmp_", tempfile(), ext))
  }

  ## Set terra temp directory once
  terra::terraOptions(tempdir = working_root, verbose = FALSE)

  ## PROCESS INPUTS ------------------------------------------------------------

  if (inherits(input, "character")) {

    if (grepl("\\.shp$", input, ignore.case = TRUE)) {
      output <- terra::vect(input)

    } else if (grepl("\\.(tif|tiff)$", input, ignore.case = TRUE)) {
      output <- terra::rast(input)

    } else {
      stop(paste0(input_name, " must end in .shp, .tif, or .tiff"))
    }

  } else if (inherits(input, c("RasterLayer", "PackedSpatRaster"))) {
    output <- terra::rast(input)

  } else if (inherits(input, c("PackedSpatVector", "sf", "sfc"))) {
    output <- terra::vect(input)

  } else if (inherits(input, c("SpatRaster", "SpatVector"))) {
    output <- input

  } else {
    stop("'input' is not a supported spatial type")
  }

  if (is.na(terra::crs(output)) || is.null(terra::crs(output))) {
    stop("'output' crs() is NULL or NA. Apply a CRS before continuing.")
  }

  orig_output <- output

  ## Subset inputs by variable names -------------------------------------------
  if (is.null(input_variable_names)) {
    input_variable_names <- names(output)
  }

  if (any(!input_variable_names %in% names(output))) {
    stop("Some 'input_variable_names' not present in input")
  }

  if (inherits(output, "SpatVector")) {
    output <- output[, input_variable_names, drop = FALSE]
  }

  if (inherits(output, "SpatRaster")) {
    output <- terra::subset(output, input_variable_names)
  }

  ## Align to align_to (if provided) -------------------------------------------
  if (!is.null(align_to)) {

    # Reuse SAME working directory for recursive calls
    align_to <- process_input(
      input = align_to,
      working_dir = working_root
    )

    if (is.na(terra::crs(align_to)) || is.null(terra::crs(align_to))) {
      stop("'align_to' crs() is NULL or NA. Apply a CRS before continuing.")
    }

    ### If align_to is a raster ------------------------------------------------
    if (inherits(align_to, "SpatRaster")) {

      ### Align vector to raster -----------------------------------------------
      if (inherits(output, "SpatVector")) {

        output <- terra::project(output, align_to)

        #### Align categorical vector to raster --------------------------------
        if (resample_type == "near") {

          # Split categorical data to multiple rasters
          output_split <- lapply(
            stats::setNames(input_variable_names, input_variable_names),
            function(varname) {
              sf_obj <- sf::st_as_sf(output)[, varname, drop = FALSE]

              splitted <- split(sf_obj, sf_obj[[varname]])
              vects <- lapply(splitted, terra::vect)

              lapply(vects, function(vx) {
                out_file <- tmp_file(".tif")
                r <- terra::rasterize(
                  x = vx, y = align_to, field = "",
                  overwrite = TRUE, ...
                )
                terra::writeRaster(r, out_file,
                                   overwrite = TRUE,
                                   gdal = "COMPRESS=NONE")
                rr <- terra::rast(out_file)
                names(rr) <- paste0(
                  varname, "_",
                  unique(vx[[varname]])[1]
                )
                rr
              })
            }
          )

          output <- terra::rast(unlist(output_split, recursive = FALSE))

        } else {

          ### Align numeric vector to raster -----------------------------------
          rasters <- lapply(input_variable_names, function(varname) {
            out_file <- tmp_file(".tif")
            r <- terra::rasterize(
              x = output,
              y = align_to,
              field = varname,
              overwrite = TRUE,
              ...
            )
            terra::writeRaster(r, out_file,
                               overwrite = TRUE,
                               gdal = "COMPRESS=NONE")
            rr <- terra::rast(out_file)
            names(rr) <- varname
            rr
          })

          output <- terra::rast(rasters)
        }
      }

      ### Align raster to raster -----------------------------------------------
      if (inherits(output, "SpatRaster")) {

        need_reproj <- !terra::compareGeom(
          x = output, y = align_to,
          lyrs = FALSE, crs = TRUE, ext = TRUE,
          rowcol = TRUE, res = TRUE,
          warncrs = FALSE, stopOnError = FALSE, messages = FALSE
        )

        if (need_reproj) {
          output <- terra::project(
            x = output,
            y = align_to,
            field = input_name,
            method = resample_type,
            overwrite = TRUE,
            ...
          )
          terra::varnames(output) <- input_name
        }
      }
    }

    ### If align_to is a vector ------------------------------------------------
    if (inherits(align_to, "SpatVector")) {

      # Raster → vector geometry conversion
      if (inherits(output, "SpatRaster")) {
        if (terra::is.polygons(align_to)) output <- terra::as.polygons(output, ...)
        if (terra::is.points(align_to))   output <- terra::as.points(output, ...)
        if (terra::is.lines(align_to))    output <- terra::as.lines(output, ...)
      }

      # Vector → vector projection
      if (inherits(output, "SpatVector")) {
        output <- terra::project(output, align_to, ...)
      }
    }
  }

  ## Clip and mask clip_region prior to heavier processing (if provided) -------
  if (!is.null(clip_region)) {

    cr <- process_input(
      input = clip_region,
      align_to = align_to,
      working_dir = working_root
    )

    if (inherits(output, "SpatVector")) {
      output <- terra::crop(output, cr)
    }

    if (inherits(output, "SpatRaster")) {
      output <- terra::crop(
        x = output, y = cr, snap = snap,
        mask = TRUE, overwrite = TRUE
      )
      output <- terra::mask(output, cr, overwrite = TRUE)
    }
  }

  ## Split categorical rasters for resample_type = "near" ----------------------
  if (inherits(output, "SpatRaster") && resample_type == "near") {

    bricks <- terra::split(output, names(output))

    bricks <- lapply(names(bricks), function(layer_name) {
      lyr <- bricks[[layer_name]]
      uvals <- unique(terra::values(lyr))
      uvals <- uvals[!is.na(uvals) & !is.nan(uvals) & !is.infinite(uvals)]

      if (length(uvals) > 1) {
        outs <- lapply(uvals, function(v) {
          vals <- terra::values(lyr)
          vals[vals != v] <- NA
          vals[vals == v] <- 1
          r <- terra::setValues(lyr, vals)
          names(r) <- paste0(layer_name, "_", v)
          r
        })
        terra::rast(outs)
      } else {
        lyr
      }
    })

    output <- terra::rast(bricks)
  }

  ## PREPARE OUTPUTS -----------------------------------------------------------

  ## If persisting final output to disk rater than in-memory -------------------

  if (isTRUE(persist_final)) {

    if (inherits(output, "SpatRaster")) {

      # If we started as a raster and geometry is identical, skip saving.
      # Otherwise, write a final file so the object is backed by a fresh TIF.
      need_save <- TRUE
      if (!inherits(orig_output, "SpatVector")) {
        same_geom <- terra::compareGeom(
          x = output, y = orig_output,
          lyrs = FALSE,   # layers not relevant for file-backing intent
          crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
          warncrs = FALSE, stopOnError = FALSE, messages = FALSE
        )
        need_save <- !isTRUE(same_geom)
      }

      if (need_save) {
        final_tif <- tmp_file(".tif")
        terra::writeRaster(output, final_tif, overwrite = TRUE, gdal = "COMPRESS=NONE")
        # Reopen to ensure the returned SpatRaster is file-backed by final_tif
        output <- terra::rast(final_tif)
      }
    }

    if (inherits(output, "SpatVector")) {
      final_shp <- tmp_file(".shp")
      terra::writeVector(output, final_shp, overwrite = TRUE)
      # Reopen to return a file-backed SpatVector
      output <- terra::vect(final_shp)
    }
  }

  ## Final output --------------------------------------------------------------
  return(output)
}
