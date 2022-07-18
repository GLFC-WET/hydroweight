#' Generate inverse distance-weighted attributes
#'
#' \code{hydroweight::hydroweight_attributes()} calculates distance-weighted
#' attributes using distance-weighted rasters generated in \code{hydroweight::hydroweight()},
#' an attribute layer (\code{loi}, e.g., land use polygon/raster), and a region of interest
#' (\code{roi}, e.g., a catchment polygon). The function outputs an attribute
#' summary table or a list that includes the summary table and layers used for calculation.
#' Summary statistics are calculated as in Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function only produces one instance of the \code{loi} x \code{distance_weights}
#' summary statistics (i.e., one \code{loi}, one \code{roi}, and one set of \code{distance_weights}).
#' See https://github.com/bkielstr/hydroweight for workflows.
#'
#' Spatial layers are aligned to \code{distance_weights} (i.e., identical coordinate reference systems - CRS).
#'
#' @param loi character (with extension, e.g., "lu.shp"), \code{sf} or \code{RasterLayer}. Layer of interest (e.g., land use layer).
#' @param loi_attr_col character. A name that will precede the calculated attributes (e.g., loi_mean, loi_median etc.)
#' @param loi_columns character. The column names over which to summarize the attributes.
#' @param loi_numeric logical. If \code{TRUE}, the `loi_columns` being summarized are numeric. If \code{FALSE}, the `loi_columns` being summarized are categorical.
#' @param loi_numeric_stats character. One or more of c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"). Those without distwtd_ are simple "lumped" statistics.
#' @param roi character (with extension, e.g., "roi.shp"), \code{sf} or \code{RasterLayer}. Region of interest (e.g., catchment boundary). Everything within this region will be used to calculate attributes.
#' @param roi_uid character. Unique identifier value for the roi.
#' @param roi_uid_col character. Column name that will be assigned to the roi_uid.
#' @param distance_weights character (with extension, e.g., "idw.rds"), \code{list}. The distance-weighted rasters output from \code{hydroweight}.
#' @param remove_region \code{sf} or \code{RasterLayer}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param return_products logical. If \code{TRUE}, a list containing attribute summary table, the \code{roi-} and \code{remove_region}-masked layer (i.e., all cells contributing to attribute calculation), and \code{distance_weights} raster. If \code{FALSE}, attribute summary table only.
#' @return If \code{return_products = TRUE}, a list containing 1) attribute summary table, and 2) a list of return_products of \code{length(distance_weights)} where each list element contains a list of 2 sub-elements: 1) \code{roi-} and \code{remove_region}-masked \code{loi} (i.e., all cells contributing to attribute calculation), and 2) the \code{roi-} and \code{remove_region}-masked \code{distance_weights} raster. If \code{return_products = FALSE}, attribute summary table only.
#' @export

hydroweight_attributes <- function(loi,
                                   loi_attr_col = NULL,
                                   loi_columns = NULL,
                                   loi_numeric = NULL,
                                   loi_numeric_stats = NULL,
                                   roi,
                                   roi_uid,
                                   roi_uid_col,
                                   distance_weights,
                                   remove_region=NULL,
                                   return_products = TRUE) {

  hydroweight_dir<-tempdir()
  if (!dir.exists(hydroweight_dir)) dir.create(hydroweight_dir)

  if (inherits(loi,"character")){
    if (grepl("\\.shp$",loi)) {
      loi<-terra::vect(loi)
    } else if (grepl("\\.tif$|\\.tiff",loi)) {
      loi<-terra::rast(loi)
    } else {
      stop("'loi' must be spcified as vector or raster layer, or a character string ending in .shp or .tif")
    }
  }
  if (inherits(loi,"RasterLayer")){
    loi<-terra::rast(loi)
  }
  if (inherits(loi,"PackedSpatRaster")){
    loi<-terra::rast(loi)
  }
  if (inherits(loi,"PackedSpatVector")){
    loi<-terra::vect(loi)
  }
  if (inherits(loi,"sf")){
    loi<-terra::vect(loi)
  }

  if (inherits(roi,"character")){
    if (grepl("\\.shp$",roi)) {
      roi<-terra::vect(roi)
    } else if (grepl("\\.tif$|\\.tiff",roi)) {
      roi<-terra::rast(roi)
    } else {
      stop("'roi' must be spcified as vector or raster layer, or a character string ending in .shp or .tif")
    }
  }
  if (inherits(roi,"RasterLayer")){
    roi<-terra::rast(roi)
  }
  if (inherits(roi,"PackedSpatRaster")){
    roi<-terra::rast(roi)
  }
  if (inherits(roi,"PackedSpatVector")){
    roi<-terra::vect(roi)
  }
  if (inherits(roi,"sf")){
    roi<-terra::vect(roi)
  }

  if (!roi_uid_col %in% names(roi)) {
    stop("`roi_uid_col` must be in `roi`")
  }

  if (inherits(remove_region,"character")){
    if (grepl("\\.shp$",remove_region)) {
      remove_region<-terra::vect(remove_region)
    } else if (grepl("\\.tif$|\\.tiff",remove_region)) {
      remove_region<-terra::rast(remove_region)
    }  else {
      stop("'remove_region' must be spcified as vector or raster layer, or a character string ending in .shp or .tif")
    }
  }
  if (inherits(remove_region,"RasterLayer")){
    remove_region<-terra::rast(remove_region)
  }
  if (inherits(remove_region,"PackedSpatRaster")){
    remove_region<-terra::rast(remove_region)
  }
  if (inherits(remove_region,"PackedSpatVector")){
    remove_region<-terra::vect(remove_region)
  }
  if (inherits(remove_region,"sf")){
    remove_region<-terra::vect(remove_region)
  }

  if (inherits(distance_weights,"character")){
    if (grepl("\\.rds$",distance_weights)) {
      distance_weights<-readRDS(distance_weights)
    } else {
      stop("'distance_weights' must be specified as .rds file")
    }
  }

  distance_weights<-lapply(distance_weights,function(x) if (inherits(x,"PackedSpatRaster")) terra::rast(x) else x)

  ## Set resampling based on loi_numeric
  if (loi_numeric == TRUE) {
    loi_resample <- "bilinear"
  }

  if (loi_numeric == FALSE) {
    loi_resample <- "ngb"
  }

  if (inherits(loi, ("SpatRaster"))) {
    message("\n Reprojecting `loi` to match `distance_weights` using method `ngb/bilinear` if loi_categorical == `TRUE/FALSE`")

    loi_r <- terra::project(
      x = loi,
      y = distance_weights[[1]],
      method = loi_resample,
      overwrite = TRUE,
      filename=file.path(paste0(tempfile(),".tif"))
    )
  }

  ## Generate RasterBrick of polygon loi using distance_weights as template

  ## If polygon loi has a field of interest with numerical data, categories
  ## are used to produce a RasterBrick of these

  ## If polygon loi has a field of interest with categorical data, categories
  ## are used to produce a RasterBrick of these
  if (inherits(loi, "SpatVector")) {
    if (loi_numeric == TRUE) {
      loi_r <- lapply(loi_columns, function(x) {
        loi_return <- terra::rasterize(x=loi,
                                       y= distance_weights[[1]],
                                       field = x
        )
      })
      names(loi_r) <- loi_columns
      loi_r <- terra::rast(loi_r)
    }

    if (loi_numeric == FALSE) {
      loi_r <- lapply(loi_columns, function(x) {
        brick_ret <- terra::rasterize(x=loi,
                                      y= distance_weights[[1]],
                                      field = x
        )
        # names(brick_ret) <- paste0(x, names(brick_ret))
        brick_ret
      })
      loi_r <- terra::rast(loi_r)

      ## Account for fasterizing NAs
      # names(loi_r) <- gsub("NA\\.", "XNA", names(loi_r))
    }
  }

  ## Generate roi
  if (inherits(roi, "SpatVector")) {
    roi_r <- terra::rasterize(x=roi, y= distance_weights[[1]])
  } else {
    message("\n Reprojecting `roi` to match `distance_weights` attributes using method `loi_resample`")

    roi_r <- terra::project(
      x = roi,
      y = distance_weights[[1]],
      method = loi_resample,
      overwrite = TRUE,
      filename=file.path(paste0(tempfile(),".tif"))
    )
  }

  ## Mask out remove_region if present
  if (!is.null(remove_region)) {
    if (inherits(remove_region, "SpatVector")) {
      remove_region_r <- terra::rasterize(x=remove_region, y = distance_weights[[1]])
    } else {
      remove_region_r <- remove_region
    }

    roi_r <- terra::mask(x=roi_r, mask = remove_region_r, inverse = TRUE)
  }

  distance_weights_attributes <- lapply(1:length(distance_weights), function(ii) {
    #browser()
    # This loop could definitely be optimized better

    ## For raster data with loi_numeric numbers
    if (loi_numeric == TRUE) {

      ## Mask out roi
      loi_r_mask <- terra::mask(loi_r, roi_r)
      distance_weights_mask <- terra::mask(distance_weights[[ii]], roi_r)

      ## Make sure NA across loi_r_mask is replicated across distance_weights_mask
      distance_weights_mask <- distance_weights_mask * !is.na(loi_r_mask)

      ##
      loi_dist <- loi_r_mask * distance_weights_mask
      names(loi_dist) <- names(loi_r_mask)

      ## Weighted mean
      loi_distwtd_mean <- terra::global(loi_dist, fun = "sum", na.rm = TRUE) / terra::global(distance_weights_mask, fun = "sum", na.rm = T)
      colnames(loi_distwtd_mean)<-"distwtd_mean"

      ## Weighted standard deviation
      ## https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
      term1 <- terra::global((distance_weights_mask * (loi_r_mask - unlist(loi_distwtd_mean))^2), fun = "sum", na.rm = TRUE)
      M <- terra::global(distance_weights_mask != 0, fun="sum", na.rm = T)
      term2 <- ((M - 1) / M) * terra::global(distance_weights_mask, fun = "sum", na.rm = TRUE)

      loi_distwtd_sd <- sqrt(term1 / term2)
      colnames(loi_distwtd_sd)<-"distwtd_sd"

      ## Non-weighted statistics
      loi_mean <- terra::global(loi_r_mask, fun = "mean", na.rm = T)
      colnames(loi_mean)<-"mean"
      loi_sd <- terra::global(loi_r_mask, fun = "sd", na.rm = T)
      colnames(loi_sd)<-"sd"
      loi_median <- terra::global(loi_r_mask, fun = function(x) median(x,na.rm=T))
      colnames(loi_median)<-"median"
      loi_min <- terra::global(loi_r_mask, fun = "min", na.rm = T)
      colnames(loi_min)<-"min"
      loi_max <- terra::global(loi_r_mask, fun = "max", na.rm = T)
      colnames(loi_max)<-"max"
      loi_sum <- terra::global(loi_r_mask, fun = "sum", na.rm = T)
      colnames(loi_sum)<-"sum"
      loi_cell_count <- !is.na(loi_r_mask)
      loi_cell_count <- terra::global(loi_cell_count, fun="sum", na.rm = T)
      colnames(loi_cell_count)<-"cell_count"
      roi_cell_count <- terra::global(roi_r, fun="sum", na.rm = T)
      loi_NA_cell_count <- unlist(roi_cell_count) - unlist(loi_cell_count)
      loi_NA_cell_count<-data.frame(matrix(loi_NA_cell_count,nrow=dim(loi_mean)[1]))
      colnames(loi_NA_cell_count)<-"NA_cell_count"

      loi_stats <- data.frame(
        loi_distwtd_mean, loi_distwtd_sd, loi_mean, loi_sd,
        loi_median, loi_min, loi_max,
        loi_sum, loi_cell_count, loi_NA_cell_count
      )

      ## If single RasterLayer
      if (terra::nlyr(loi_dist) == 1) {
        loi_stats <- as.data.frame(loi_stats)
        colnames(loi_stats) <- gsub("loi_", "", colnames(loi_stats))
        colnames(loi_stats) <- paste(loi_attr_col, colnames(loi_stats), sep = "_")
      }

      ## If RasterBrick
      if (terra::nlyr(loi_dist) > 1) {
        loi_stats$cats <- rownames(loi_stats)
        loi_stats_w <- tidyr::pivot_wider(loi_stats, values_from = c(any_of(colnames(loi_stats))), names_from = "cats")
        colnames(loi_stats_w) <- gsub("loi_", "", colnames(loi_stats_w))

        vars_strings <- stringr::str_extract(colnames(loi_stats_w), loi_columns)
        stats_strings <- stringr::str_replace(colnames(loi_stats_w), paste0("_", loi_columns), "")
        colnames_strings <- paste(vars_strings, stats_strings, sep = "_")

        loi_stats_w <- as.data.frame(loi_stats_w)
        colnames(loi_stats_w) <- paste(loi_attr_col, colnames_strings, sep = "_")
        loi_stats_w <- as.data.frame(loi_stats_w)

        loi_stats <- loi_stats_w
      }
    }

    ## For raster data with categorical numbers
    if (loi_numeric == FALSE) {

      ## Construct brick if "RasterLayer"
      if (inherits(loi_r, "SpatRaster")) {

        if (nlyr(loi_r)>1) {
          brick_list<-terra::split(loi_r,names(loi_r))
          names(brick_list)<-names(loi_r)
        } else {
          brick_list<-loi_r
        }

        brick_list<-lapply(names(brick_list),function(y){
          uv <- unlist(terra::unique(brick_list[[y]]))

          out<-lapply(uv,function(x) {
            repl<-values(brick_list[[y]])
            repl[repl!=x]<-NA
            repl[!is.na(repl)]<-1
            out<-terra::setValues(brick_list[[y]],repl)
            names(out)<-paste0(y,"_",x)
            return(out)
          })

          nms<-sapply(out,names)
          out <- terra::rast(out)
          names(out) <- nms
          return(out)
        })

        loi_r <- terra::rast(brick_list)

        #names(loi_r) <- uv


      }

      loi_r_mask <- terra::mask(loi_r, roi_r)
      distance_weights_mask <- terra::mask(distance_weights[[ii]], roi_r)

      loi_dist <- loi_r_mask * distance_weights_mask
      names(loi_dist) <- names(loi_r_mask)

      loi_pct_distwtd <- unlist(terra::global(loi_dist, fun = "sum", na.rm = TRUE)) /
        unlist(terra::global(distance_weights_mask, fun = "sum", na.rm = TRUE))
      names(loi_pct_distwtd) <- names(loi_r_mask)

      loi_stats <- data.frame(t(loi_pct_distwtd))
      colnames(loi_stats) <- gsub("X", "_", colnames(loi_stats))
      colnames(loi_stats) <- paste(loi_attr_col, "prop", colnames(loi_stats), sep = "_")
      colnames(loi_stats) <- gsub("__", "_", colnames(loi_stats))
    }

    loi_stats <- data.frame(data.frame(roi_uid = roi_uid), loi_stats)
    colnames(loi_stats)[1] <- roi_uid_col

    ## Reduce loi_stats frame according to loi_numeric_stats
    if (!is.null(loi_numeric_stats)) {
      col_return <- lapply(loi_numeric_stats, function(x) {
        grep(x, names(loi_stats))
      })
      col_return <- do.call("c", col_return)
      col_return <- unique(col_return)

      loi_stats <- loi_stats[, c(1, col_return)]
    }

    #loi_stats

    if (return_products == TRUE) {
      ret_list <- list(loi_stats, loi_r_mask, distance_weights_mask)
      names(ret_list) <- c("loi_statistics", "loi_Raster*_bounded", "distance_weights_bounded")
    } else {
      ret_list <- loi_stats
    }

    return(ret_list)
  })

  ## Pull out attribute_frames and adjust
  attribute_frames <- lapply(1:length(distance_weights_attributes), function(ii) {
    if (return_products == TRUE) {
      sub_frame <- distance_weights_attributes[[ii]][[1]]

      colnames(sub_frame)[-1] <- paste0(
        names(distance_weights)[[ii]],
        "_",
        colnames(sub_frame)[-1]
      )
    }

    if (return_products == FALSE) {
      sub_frame <- distance_weights_attributes[[ii]]

      colnames(sub_frame)[-1] <- paste0(
        names(distance_weights)[[ii]],
        "_",
        colnames(sub_frame)[-1]
      )
    }

    return(sub_frame)
  })
  attribute_frames <- Reduce(merge, attribute_frames)

  if (loi_numeric == TRUE) {

    ## Pull out distance-weighted columns
    distwtd_cols <- grep("distwtd", names(attribute_frames))
    distwtd_data <- attribute_frames[distwtd_cols]

    ## All lumped columns are identical data regardless of prefix. Reducing.
    lumped_cols <- which(!grepl("distwtd", colnames(attribute_frames)) & !grepl(roi_uid_col, colnames(attribute_frames)))
    lumped_data <- attribute_frames[lumped_cols]
    lumped_cols_reduce <- grep(names(distance_weights)[1], colnames(lumped_data))
    lumped_data <- lumped_data[lumped_cols_reduce]
    colnames(lumped_data) <- gsub(names(distance_weights)[1], "lumped", colnames(lumped_data))

    ## Bring together
    attribute_frames_return <- cbind(
      attribute_frames[, roi_uid_col],
      lumped_data,
      distwtd_data
    )
    colnames(attribute_frames_return)[1] <- roi_uid_col
  } else {
    attribute_frames_return <- attribute_frames
  }

  ## If return_products == TRUE, remove the attribute table, keep the other parts of the list, and name according to distance_weights
  if (return_products == TRUE) {
    attribute_products <- lapply(1:length(distance_weights_attributes), function(ii) {
      sub_products <- lapply(distance_weights_attributes[[ii]][-1],terra::wrap)
    })
    names(attribute_products) <- names(distance_weights)

    return_list <- list(attribute_frames_return, attribute_products)
    names(return_list) <- c("attribute_table", "return_products")
  } else {
    return_list <- attribute_frames_return
  }

  return(return_list)
}

