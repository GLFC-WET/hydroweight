#' Align and normalize spatial inputs for hydroweight workflows
#'
#' `process_input()` standardizes spatial inputs (vectors or rasters, or
#' filepaths to either) into `terra` objects and optionally:
#' (1) subsets requested variables/layers, (2) aligns to a target geometry,
#' and (3) clips/masks to a region of interest. The result is a `SpatVector`
#' or `SpatRaster` ready for downstream use in `hydroweight()` and
#' `hydroweight_attributes()`.
#'
#' @section Coercion & alignment (high‑level):
#' * Vector or raster inputs (or filepaths) are coerced to `SpatVector` or
#'   `SpatRaster` as appropriate.
#' * If `align_to` is a raster, rasters are reprojected/resampled to match its
#'   CRS, extent, resolution, and grid. Vectors are projected and, when needed,
#'   rasterized against the `align_to` grid.
#' * If `resample_type = "near"`, categorical data are handled using nearest‑neighbor
#'   semantics; if `"bilinear"`, continuous data are interpolated bilinearly.
#' * When `clip_region` is supplied, outputs are cropped and (for rasters) masked
#'   to that region, with grid placement controlled by `snap`.
#'
#' @param input May be a spatial object (`sf`, `SpatVector`, `RasterLayer`, `SpatRaster`, etc.) or a filepath
#'   to `.shp` dataset readable by `sf::st_read()` (i.e., Shapefile)
#'   or `.tif` dataset readable by `terra::rast()` (i.e., GeoTIFF).
#'
#' @param input_name Optional label (`character`) used in diagnostic messages.
#'   Helpful for tracing errors when processing multiple inputs.
#'
#' @param input_variable_names Optional `character` vector of attribute names
#'   (for vectors) or layer names (for rasters) to retain. If `NULL`, all
#'   attributes/layers are kept.
#'
#' @param align_to Optional **target** to align to: a spatial object or filepath
#'   (same accepted formats as `input`). If supplied, the returned object is
#'   projected/resampled to the CRS and grid of `align_to` (for rasters) or
#'   projected to its CRS (for vectors). If `NULL`, inputs are returned in
#'   their native geometry after optional subsetting.
#'
#' @param clip_region Optional spatial object or filepath (same accepted formats as `input`).
#'   If provided, the result is cropped to `clip_region`; rasters are also masked.
#'
#' @param resample_type Character; resampling strategy for raster alignment and/or
#'   rasterization to `align_to`. One of `"bilinear"` (typical for continuous
#'   data) or `"near"` (typical for categorical data). Default: `"bilinear"`.
#'
#' @param snap Character. Controls how `terra::crop()` aligns the raster grid to the
#'   clipping geometry (`"near"`, `"in"`, or `"out"`). The choice of snap setting
#'   controls how raster grids align during cropping, which can shift cell boundaries
#'   and influence flow‑path continuity, ultimately affecting hydrological distance
#'   and accumulation calculations. (see *Snap Argument* section below)
#'
#' @param working_dir Optional directory path used for temporary files and
#'   `terra`’s out‑of‑memory operations. If `NULL`, a temporary folder is created.
#'
#' @param persist_final Logical; if `TRUE`, the final `SpatRaster`/`SpatVector`
#'   is written to disk inside `working_dir` and the function returns a
#'   file‑backed object. Default: `FALSE`.
#'
#' @param ... Additional arguments forwarded to underlying `terra` operations
#'   (e.g., `terra::crop()`, `terra::project()`, `terra::rasterize()`,
#'   `terra::as.polygons()`).
#'
#' @return A `SpatVector` or `SpatRaster` (file‑backed if `persist_final = TRUE`),
#'   aligned and optionally subsetted and clipped as requested.
#'
#' @section Notes:
#' * `resample_type` controls resampling when aligning to raster targets:
#'   use `"bilinear"` for numeric/continuous layers and `"near"` for categorical.
#' * `snap` controls grid placement during cropping and can influence flow‑path
#'   continuity in hydrologic analyses.
#' * When aligning a vector to a raster target, categorical attributes may be
#'   rasterized into multiple binary layers (one per category) when
#'   `resample_type = "near"`.
#'
#' @seealso
#' * [hydroweight()] for generating distance‑weighted rasters
#' * [hydroweight_attributes()] for distance‑weighted attribute summaries
#'
#' @export
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
    ...) {

  ## SETUP ---------------------------------------------------------------------

  if (is.null(input)) {
    return(input)
  }

  ## Allowed choices
  valid_rt   <- c("bilinear", "near")
  valid_snap <- c("near","in","out")

  resample_type <- match.arg(resample_type, choices = valid_rt)
  snap <- match.arg(snap, choices = valid_snap, several.ok = FALSE)

  ## Use a single stable working directory per call
  working_root <- if (is.null(working_dir)) {
    tempfile("hwt_")
  } else {
    working_dir
  }
  dir.create(working_root, showWarnings = FALSE, recursive = TRUE)

  ## Add a local terraOptions guard and restore it on exit
  ## Each call to process_input() now isolates its tempdir change and never
  ## contaminates the rest of the worker, which is crucial when process_input()
  ## is called recursively or from different places that expect a consistent temp root.
  old_terra_opts <- terra::terraOptions(print = FALSE)
  on.exit({
    opts <- old_terra_opts
    opts$metadata <- NULL; opts$names <- NULL
    if (identical(opts$filetype, "")) opts$filetype <- NULL
    do.call(terra::terraOptions, opts)
  }, add = TRUE)

  ## Set terra temp directory once for the duration of this call
  terra::terraOptions(tempdir = working_root, verbose = FALSE)

  ## Helper for generating temporary file paths
  tmp_file <- function(ext = ".tif") {
    file.path(working_root, paste0("tmp_", basename(tempfile()), ext))
  }


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

  if (inherits(output, "SpatVector")) {
    if (length(input_variable_names) > 0) {
      output <- output[, input_variable_names, drop = FALSE]
    }
  }

  if (any(!input_variable_names %in% names(output))) {
    stop("Some 'input_variable_names' not present in input")
  }

  if (inherits(output, "SpatRaster")) {
    output <- terra::subset(output, input_variable_names)
  }

  ## Align to align_to (if provided) -------------------------------------------
  if (!is.null(align_to)) {

    # Reuse SAME working directory for recursive calls
    align_to <- process_input(
      input = align_to,
      input_name = "align_to",
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
                  overwrite = TRUE,
                  ...
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
        x = output, y = cr, snap = "near",
        mask = TRUE, overwrite = TRUE
      )
      output <- terra::mask(output, cr, overwrite = TRUE)
    }
  }

  ## Split categorical rasters for resample_type = "near" ----------------------
  if (inherits(output, "SpatRaster") && resample_type == "near") {

    if(terra::nlyr(output) > 1){
      brick_list<-terra::split(output,names(output))
      names(brick_list)<-names(output)
    } else {
      brick_list<-output
    }

    brick_list <- lapply(names(brick_list), function(y){

          b <- as.numeric(brick_list[[y]])
          f <- brick_list[[y]]
          (uv <- unlist(terra::unique(b)))
          (uv <- uv[!sapply(uv,is.nan)])
          (uv <- uv[!sapply(uv,is.na)])
          (uv <- uv[!sapply(uv,is.infinite)])

          if (length(uv)>1) {
            out<-lapply(uv,function(x) {
              repl <- terra::values(b)
              repl[repl!=x] <- NA
              repl[!is.na(repl)]<-1
              out<-terra::setValues(b,repl)
              if (length(uv)>1) names(out) <-paste0(y,"_",x)
              if (length(uv)==1) names(out) <- y
              return(out)
            })

            nms <-sapply(out, names)
            out <- terra::rast(out)
            names(out) <- nms



          } else {
            out<-brick_list[[y]]
          }
          return(out)
        })

        output <- terra::rast(brick_list)
      }

  #   bricks <- terra::split(output, names(output))
  #   names(bricks) <- names(output)
  #
  #   bricks <- lapply(names(bricks), function(layer_name) {
  #     lyr <- bricks[[layer_name]]
  #
  #     # Pull factor levels (if present) to retrieve human-readable labels
  #     levs_list <- if (terra::is.factor(lyr)) levels(lyr) else NULL  # base::levels -> terra method
  #     levs_tbl  <- if (!is.null(levs_list) && length(levs_list) >= 1) levs_list[[1]] else NULL
  #     if (!is.null(levs_tbl)) {
  #       id_col  <- if ("ID" %in% names(levs_tbl)) "ID" else names(levs_tbl)[1]
  #       lab_col <- setdiff(names(levs_tbl), id_col)
  #       lab_col <- if (length(lab_col) > 0) lab_col[1] else id_col
  #     }
  #
  #     # Present category codes in this (already clipped) layer
  #     f <- try(terra::freq(lyr), silent = TRUE)     # NOTE: terra::freq() has no 'useNA' arg
  #     uvals <- if (!inherits(f, "try-error") && !is.null(f) && nrow(f) > 0) {
  #       f[, "value"]
  #     } else {
  #       unique(terra::values(lyr))
  #     }
  #     uvals <- uvals[is.finite(uvals)]
  #
  #     if (length(uvals) >= 1) {
  #       # IMPORTANT: also handle the single-class case (length == 1)
  #       outs <- lapply(uvals, function(v) {
  #         # Make a binary mask (1/NA) so proportions compute correctly
  #         r <- terra::ifel(lyr == v, 1L, NA)
  #
  #         # Convert to factor, attach a one-row levels() with ID=1 and the label for this class
  #         r <- terra::as.factor(r)
  #         lab <- if (!is.null(levs_tbl)) {
  #           L <- levs_tbl[[lab_col]][levs_tbl[[id_col]] == v]
  #           if (length(L) == 1 && !is.na(L)) as.character(L) else as.character(v)
  #         } else {
  #           as.character(v)
  #         }
  #         levels(r) <- list(data.frame(ID = 1L, Class = lab, stringsAsFactors = FALSE))
  #
  #         # Name with label for readability (fallback to code if no label)
  #         safe_lab <- gsub("[^A-Za-z0-9_]+", "_", lab)
  #         names(r) <- paste0(layer_name, "_", safe_lab)
  #         r
  #       })
  #
  #       terra::rast(outs)
  #     } else {
  #       # Nothing present (all NA) — pass through unchanged
  #       lyr
  #     }
  #   })
  #
  #   output <- terra::rast(bricks)
  # }

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
