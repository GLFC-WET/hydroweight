#' Generate inverse distance-weighted attributes
#'
#' \code{hydroweight::hydroweight_attributes()} calculates distance-weighted attributes using distance-weighted rasters generated in \code{hydroweight::hydroweight()}, a region of interest
#' (\code{roi}, e.g., a catchment polygon), and an attribute layer (\code{loi}, e.g., land use polygon/raster). The function outputs the attribute
#' summary table or a list including the summary table and layers used for calculation. Summary statistics are calculated as in Peterson et al. 2011 <https://doi:10.1111/j.1365-2427.2010.02507.x>).
#' IMPORTANTLY, this function only produces one instance of the layer x weight summary statistics (i.e., one \code{roi}, one \code{loi}, and one \code{distance_weight}).
#' See example and vignette for simple workflow.
#'
#' Spatial layers should align (i.e., identical coordinate reference systems - CRS).
#'
#' @param roi \code{sf} or \code{RasterLayer}. Region of interest (e.g., catchment boundary). Everything within this region will be used to calculate attributes.
#' @param loi \code{sf} or \code{RasterLayer}. Layer of interest (e.g., land use layer).
#' @param loi_numeric logical. If \code{TRUE}, the attributes being summarized are numeric. If \code{FALSE}, the attributes being summarized are categorical.
#' @param loi_categories character. If \code{loi_numeric = FALSE}, the column names over which to summarize the attributes.
#' @param distance_weight \code{RasterLayer}. The distance-weighted raster.
#' @param remove_region \code{sf} or \code{RasterLayer}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param resample character; \code{ngb} or \code{bilinear}. If re-projection needed to distance_weight projection, should nearest neighbour (categorical: "ngb") or bilinear (numeric: "bilinear") resampling be used?
#' @param return_products logical. If \code{TRUE}, a list containing attribute summary table, the \code{roi-} and \code{remove_region}-masked layer (i.e., all cells contributing to attribute calculation), and \code{distance_weight} raster. If \code{FALSE}, attribute summary table only.
#' @param uid character. Unique identifier value for the roi
#' @param uid_col character. Column name that will be assigned to the uid.
#' @return If \code{return_products = TRUE}, a list containing attribute summary table, the \code{roi-} and \code{remove_region}-masked layer (i.e., all cells contributing to attribute calculation), and \code{distance_weight} raster. If \code{return_products = FALSE}, attribute summary table only.
#' @export

hydroweight_attributes <- function(roi = NULL,
                                   loi = NULL,
                                   loi_numeric = NULL,
                                   loi_categories = NULL,
                                   distance_weight = NULL,
                                   remove_region = NULL,
                                   resample = NULL,
                                   return_products = TRUE,
                                   uid = NULL,
                                   uid_col = NULL)
{

  if(class(loi)[1] == "RasterLayer"){

      message("Reprojecting to distance_weight extent, resolution, and origin.
      This could be very slow given our choice of algorithm.
      Consider reprojecting loi before this step.")

        loi_r <- raster::projectRaster(from = loi,
                                       to  = distance_weight,
                                       method = resample)

  }

  ## Generate RasterBrick of polygon loi using distance_weight as template

  ## If polygon loi has a field of interest with numerical data, categories
  ## are used to produce a RasterBrick of these

  ## If polygon loi has a field of interest with categorical data, categories
  ## are used to produce a RasterBrick of these

  if(class(loi)[1] =="sf"){

    if(loi_numeric == TRUE){

      loi_r <- lapply(categories, function(x){

        loi_return <- fasterize::fasterize(loi, raster = distance_weight,
                                           field = x)
      })
      names(loi_r) <- categories
      loi_r <- brick(loi_r)
    }

    if(loi_numeric == FALSE){

      loi_r <- lapply(categories, function(x){

        brick_ret <- fasterize::fasterize(loi, raster = distance_weight,
                                          by = x)
        names(brick_ret) <- paste0(x, names(brick_ret))
        brick_ret
      })
      loi_r <- brick(loi_r)
    }

  }

  ## Mask to roi
  if(class(roi)[1] == "sf"){
    roi_r <-  fasterize::fasterize(roi, raster = distance_weight)
  } else {
    roi_r <- roi
  }

  loi_r_mask <- raster::mask(loi_r, roi_r)
  distance_weight_mask <- raster::mask(distance_weight, roi_r)

  ## Mask out region
  if(!is.null(remove_region)){

    if(class(remove_region)[1] == "sf"){
      remove_region_r <-  fasterize::fasterize(remove_region, raster = distance_weight)
    } else {
      remove_region_r <- remove_region
    }

    loi_r_mask <- raster::mask(loi_r_mask, remove_region_r,
                               inverse = TRUE)
    distance_weight_mask <- raster::mask(distance_weight_mask, remove_region_r,
                                         inverse = TRUE)
  }

  ## For raster data with loi_numeric numbers
  if(loi_numeric == TRUE){

    (loi_dist <- loi_r_mask * distance_weight_mask)
    names(loi_dist) <- names(loi_r_mask)

    loi_dist_ret <- lapply(names(loi_dist), function(x){

      (loi_mean <- mean(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_median <- median(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_sd <- sd(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_sum <- sum(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_min <- min(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_max <- max(loi_r_mask[[x]]@data@values, na.rm=T))
      (loi_pixel_count <- sum(!is.na(loi_r_mask[[x]]@data@values)))

      ## Weighted mean
      (loi_mean_distwtd <- loi_sum / sum(distance_weight_mask@data@values, na.rm=T))

      ## Weighted standard deviation
      loi_sd_distwtd <- "NotCalculated"
      ## Not sure how to calculate yet

      loi_stats <- data.frame(loi_mean, loi_median, loi_sd, loi_sum,
                              loi_min, loi_max,
                              loi_mean_distwtd, loi_sd_distwtd,
                              loi_pixel_count)
      colnames(loi_stats) <- gsub("loi", x, colnames(loi_stats))
      return(loi_stats)
    })

    loi_stats <- do.call(dplyr::bind_cols, loi_dist_ret)

  }

  ## For raster data with categorical numbers
  if(loi_numeric == FALSE){

    ## Construct brick if "RasterLayer"
    if(class(loi_r_mask) == "RasterLayer"){

      (uv <- unique(loi_r_mask@data@values))
      (uv <- uv[!is.na(uv)])

      brick_list <- lapply(uv, function(x){

        loi_r_mask_ret <- loi_r_mask
        loi_r_mask_ret[loi_r_mask_ret@data@values == x] <- 9999
        loi_r_mask_ret[loi_r_mask_ret@data@values != 9999] <- NA
        loi_r_mask_ret[loi_r_mask_ret@data@values == 9999] <- 1

        return(loi_r_mask_ret)
      })

      loi_r_mask <- brick(brick_list)
      names(loi_r_mask) <- uv
    }

    loi_dist <- loi_r_mask * distance_weight_mask
    names(loi_dist) <- names(loi_r_mask)

    (loi_pct_distwtd <- cellStats(loi_dist, stat = "sum") /
        cellStats(distance_weight_mask, stat = "sum"))
    (names(loi_pct_distwtd) <- names(loi_r_mask))

    (loi_stats <- data.frame(t(loi_pct_distwtd)))

  }

  loi_stats$uid <- uid

  if(return_products == TRUE){

    ret_list <- list(loi_stats, loi_r_mask, distance_weight_mask)
    names(ret_list) <- c("loi_statistics", "loi_Raster*_bounded", "distance_weight_bounded")

  } else {

    ret_list <- loi_stats

  }

  return(ret_list)

}
