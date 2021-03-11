#' Generate inverse distance-weighted attributes
#'
#' \code{hydroweight::hydroweight_attributes()} calculates distance-weighted
#' attributes using distance-weighted rasters generated in \code{hydroweight::hydroweight()},
#' an attribute layer (\code{loi}, e.g., land use polygon/raster), and a region of interest
#' (\code{roi}, e.g., a catchment polygon). The function outputs an attribute
#' summary table or a list that includes the summary table and layers used for calculation.
#' Summary statistics are calculated as in Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function only produces one instance of the \code{loi} x \code{distance_weight}
#' summary statistics (i.e., one \code{loi}, one \code{roi}, and one \code{distance_weight}).
#' See vignette for workflows.
#'
#' Spatial layers are aligned to \code{distance_weight} (i.e., identical coordinate reference systems - CRS).
#'
#' @param loi \code{sf} or \code{RasterLayer}. Layer of interest (e.g., land use layer).
#' @param loi_attr_col character. A name that will precede the attributes (e.g., loi_mean, loi_median etc.)
#' @param loi_categories character. If \code{loi_numeric = FALSE}, the column names over which to summarize the attributes.
#' @param loi_numeric logical. If \code{TRUE}, the attributes being summarized are numeric. If \code{FALSE}, the attributes being summarized are categorical.
#' @param loi_numeric_stats character. One or more of c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "pixel_count"). Those without distwtd_ are simple "lumped" statistics.
#' @param roi \code{sf} or \code{RasterLayer}. Region of interest (e.g., catchment boundary). Everything within this region will be used to calculate attributes.
#' @param roi_uid character. Unique identifier value for the roi.
#' @param roi_uid_col character. Column name that will be assigned to the roi_uid.
#' @param distance_weight \code{RasterLayer}. The distance-weighted raster.
#' @param remove_region \code{sf} or \code{RasterLayer}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param return_products logical. If \code{TRUE}, a list containing attribute summary table, the \code{roi-} and \code{remove_region}-masked layer (i.e., all cells contributing to attribute calculation), and \code{distance_weight} raster. If \code{FALSE}, attribute summary table only.
#' @return If \code{return_products = TRUE}, a list containing 1) attribute summary table, 2) the \code{roi-} and \code{remove_region}-masked \code{loi} (i.e., all cells contributing to attribute calculation), and 3) the \code{roi-} and \code{remove_region}-masked \code{distance_weight} raster. If \code{return_products = FALSE}, attribute summary table only.
#' @export

hydroweight_attributes <- function(loi = NULL,
                                   loi_attr_col = NULL,
                                   loi_categories = NULL,
                                   loi_numeric = NULL,
                                   loi_numeric_stats= NULL,
                                   roi = NULL,
                                   roi_uid = NULL,
                                   roi_uid_col = NULL,
                                   distance_weight = NULL,
                                   remove_region = NULL,
                                   return_products = TRUE) {

  ## Set resampling based on loi_numeric
  if (loi_numeric == TRUE){

    loi_resample <- "bilinear"

  }

  if (loi_numeric == FALSE){

    loi_resample <- "ngb"

  }

  if (class(loi)[1] == "RasterLayer") {
    message("\n Reprojecting `loi` to match `distance_weight` using method `ngb/bilinear` if loi_categorical == `TRUE/FALSE`")

    loi_r <- raster::projectRaster(
      from = loi,
      to = distance_weight,
      method = loi_resample
    )
  }

  ## Generate RasterBrick of polygon loi using distance_weight as template

  ## If polygon loi has a field of interest with numerical data, categories
  ## are used to produce a RasterBrick of these

  ## If polygon loi has a field of interest with categorical data, categories
  ## are used to produce a RasterBrick of these

  if (class(loi)[1] == "sf") {
    if (loi_numeric == TRUE) {
      loi_r <- lapply(loi_categories, function(x) {
        loi_return <- fasterize::fasterize(loi,
          raster = distance_weight,
          field = x
        )
      })
      names(loi_r) <- loi_categories
      loi_r <- raster::brick(loi_r)
    }

    if (loi_numeric == FALSE) {
      loi_r <- lapply(loi_categories, function(x) {
        brick_ret <- fasterize::fasterize(loi,
          raster = distance_weight,
          by = x
        )
        names(brick_ret) <- paste0(x, names(brick_ret))
        brick_ret
      })
      loi_r <- raster::brick(loi_r)
    }
  }

  ## Mask to roi
  if (class(roi)[1] == "sf") {
    roi_r <- fasterize::fasterize(roi, raster = distance_weight)
  } else {
    message("\n Reprojecting `roi` to match `distance_weight` attributes using method `loi_resample`")

    roi_r <- raster::projectRaster(
      from = roi,
      to = distance_weight,
      method = loi_resample
    )
  }

  loi_r_mask <- raster::mask(loi_r, roi_r)
  distance_weight_mask <- raster::mask(distance_weight, roi_r)

  ## Mask out remove_region
  if (!is.null(remove_region)) {
    if (class(remove_region)[1] == "sf") {
      remove_region_r <- fasterize::fasterize(remove_region, raster = distance_weight)
    } else {
      remove_region_r <- remove_region
    }

    loi_r_mask <- raster::mask(loi_r_mask, remove_region_r,
      inverse = TRUE
    )
    distance_weight_mask <- raster::mask(distance_weight_mask, remove_region_r,
      inverse = TRUE
    )
  }

  ## For raster data with loi_numeric numbers
  if (loi_numeric == TRUE) {
    (loi_dist <- loi_r_mask * distance_weight_mask)
    names(loi_dist) <- names(loi_r_mask)

    ## Weighted mean
    (loi_distwtd_mean <- raster::cellStats(loi_dist, stat = "sum", na.rm = TRUE) / raster::cellStats(distance_weight_mask, stat = "sum", na.rm = T))

    ## Weighted standard deviation
    ## https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
    (term1 <- raster::cellStats((distance_weight_mask * (loi_r_mask - loi_distwtd_mean)^2), stat = "sum", na.rm = TRUE))
    (M <- raster::cellStats(distance_weight_mask != 0, "sum", na.rm = T))
    (term2 <- ((M - 1) / M) * raster::cellStats(distance_weight_mask, stat = "sum", na.rm = TRUE))

    (loi_distwtd_sd <- sqrt(term1 / term2))

    ## Non-weighted statistics
    (loi_mean <- raster::cellStats(loi_r_mask, stat = "mean", na.rm = T))
    (loi_sd <- raster::cellStats(loi_r_mask, stat = "sd", na.rm = T))
    (loi_median <- stats::median(loi_r_mask@data@values, na.rm = T))
    (loi_min <- raster::cellStats(loi_r_mask, stat = "min", na.rm = T))
    (loi_max <- raster::cellStats(loi_r_mask, stat = "max", na.rm = T))
    (loi_sum <- raster::cellStats(loi_r_mask, stat = "sum", na.rm = T))
    (loi_pixel_is_na <- !is.na(loi_r_mask))
    (loi_pixel_count <- raster::cellStats(loi_pixel_is_na, "sum", na.rm = T))

    loi_stats <- data.frame(
      loi_distwtd_mean, loi_distwtd_sd, loi_mean, loi_sd,
      loi_median, loi_min, loi_max,
      loi_sum, loi_pixel_count
    )

    ## If single RasterLayer
    if (raster::nlayers(loi_dist) == 1) {
      loi_stats <- as.data.frame(loi_stats)
      colnames(loi_stats) <- gsub("loi_", "", colnames(loi_stats))
      colnames(loi_stats) <- paste(loi_attr_col, colnames(loi_stats), sep = "_")
    }

    ## If RasterBrick
    if (raster::nlayers(loi_dist) > 1) {
      loi_stats$cats <- rownames(loi_stats)
      loi_stats_w <- tidyr::pivot_wider(loi_stats, values_from = c(
        loi_distwtd_mean, loi_distwtd_sd, loi_mean, loi_sd,
        loi_median, loi_min, loi_max,
        loi_sum, loi_pixel_count
      ), names_from = "cats")
      colnames(loi_stats_w) <- gsub("loi_", "", colnames(loi_stats_w))

      (vars_strings <- stringr::str_extract(colnames(loi_stats_w), loi_categories))
      (stats_strings <- stringr::str_replace(colnames(loi_stats_w), paste0("_", loi_categories), ""))
      (colnames_strings <- paste(vars_strings, stats_strings, sep = "_"))

      loi_stats_w <- as.data.frame(loi_stats_w)
      colnames(loi_stats_w) <- paste(loi_attr_col, colnames_strings, sep = "_")
      loi_stats_w <- as.data.frame(loi_stats_w)

      loi_stats <- loi_stats_w
    }
  }

  ## For raster data with categorical numbers
  if (loi_numeric == FALSE) {

    ## Construct brick if "RasterLayer"
    if (class(loi_r_mask) == "RasterLayer") {
      (uv <- unique(loi_r_mask@data@values))
      (uv <- uv[!is.na(uv)])

      brick_list <- lapply(uv, function(x) {
        loi_r_mask_ret <- loi_r_mask
        loi_r_mask_ret[loi_r_mask_ret@data@values == x] <- 9999
        loi_r_mask_ret[loi_r_mask_ret@data@values != 9999] <- NA
        loi_r_mask_ret[loi_r_mask_ret@data@values == 9999] <- 1

        return(loi_r_mask_ret)
      })

      loi_r_mask <- raster::brick(brick_list)
      names(loi_r_mask) <- uv
    }

    loi_dist <- loi_r_mask * distance_weight_mask
    names(loi_dist) <- names(loi_r_mask)

    (loi_pct_distwtd <- raster::cellStats(loi_dist, stat = "sum", na.rm = TRUE) /
      raster::cellStats(distance_weight_mask, stat = "sum", na.rm = TRUE))
    (names(loi_pct_distwtd) <- names(loi_r_mask))

    (loi_stats <- data.frame(t(loi_pct_distwtd)))
    (colnames(loi_stats) <- gsub("X", "_", colnames(loi_stats)))
    (colnames(loi_stats) <- paste(loi_attr_col, "prop", colnames(loi_stats), sep = "_"))
    (colnames(loi_stats) <- gsub("__", "_", colnames(loi_stats)))
  }

  loi_stats <- data.frame(data.frame(roi_uid = roi_uid), loi_stats)
  colnames(loi_stats)[1] <- roi_uid_col
  colnames(loi_stats) <- gsub("inv_", "", colnames(loi_stats))

  ## Reduce loi_stats frame according to loi_numeric_stats
  if(!is.null(loi_numeric_stats)){

    col_return <- lapply(loi_numeric_stats, function(x){grep(x, names(loi_stats))})
    col_return <- do.call("c", col_return)
    col_return <- unique(col_return)

    loi_stats <- loi_stats[,c(1, col_return)]
  }

  if (return_products == TRUE) {
    ret_list <- list(loi_stats, loi_r_mask, distance_weight_mask)
    names(ret_list) <- c("loi_statistics", "loi_Raster*_bounded", "distance_weight_bounded")
  } else {
    ret_list <- loi_stats
  }

  return(ret_list)
}
