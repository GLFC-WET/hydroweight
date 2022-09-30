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
#' @param loi character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Layer of interest (e.g., land use layer).
#' @param loi_columns character. The column names over which to summarize the attributes.
#' @param loi_numeric logical. If \code{TRUE}, the `loi_columns` being summarized are numeric. If \code{FALSE}, the `loi_columns` being summarized are categorical.
#' @param loi_numeric_stats character. One or more of c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"). Those without distwtd_ are simple "lumped" statistics.
#' @param roi character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/roi.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Region of interest (e.g., catchment boundary). Everything within this region will be used to calculate attributes.
#' @param roi_uid character. Unique identifier value for the roi.
#' @param roi_uid_col character. Column name that will be assigned to the roi_uid; Will be included in the attribute summary table.
#' @param distance_weights character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/idw.rds", or "idw.zip"), \code{list}. The distance-weighted rasters output from \code{hydroweight}.
#' @param remove_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Regions to remove when summarizing the attributes (e.g., remove lake from catchment)
#' @param return_products logical. If \code{TRUE}, a list containing attribute summary table, the \code{roi-} and \code{remove_region}-masked layer (i.e., all cells contributing to attribute calculation), and \code{distance_weights} raster. If \code{FALSE}, attribute summary table only.
#'
#' @return If \code{return_products = TRUE}, a list containing 1) attribute summary table, and 2) a list of return_products of \code{length(distance_weights)} where each list element contains a list of 2 sub-elements: 1) \code{roi-} and \code{remove_region}-masked \code{loi} (i.e., all cells contributing to attribute calculation), and 2) the \code{roi-} and \code{remove_region}-masked \code{distance_weights} raster. If \code{return_products = FALSE}, attribute summary table only.
#' @export

hydroweight_attributes <- function(loi,
                                   loi_columns = NULL,
                                   loi_numeric,
                                   loi_numeric_stats = c("distwtd_mean", "distwtd_sd", "mean", "sd", "median", "min", "max", "sum", "cell_count"),
                                   roi=NULL,
                                   roi_uid,
                                   roi_uid_col="ID",
                                   distance_weights,
                                   remove_region=NULL,
                                   return_products = TRUE) {
  require(terra)

  loi_numeric_stats<-match.arg(loi_numeric_stats,several.ok = T)

  if (loi_numeric == TRUE) {
    loi_resample <- "bilinear"
  }

  if (loi_numeric == FALSE) {
    loi_resample <- "near"
  }

  # Read in distance weights ------------------------------------------------
  own_tempdir<-gsub("file","",tempfile())
  if (!dir.exists(own_tempdir)) dir.create(own_tempdir)
  terra::terraOptions(tempdir = own_tempdir, verbose=F)

  if (inherits(distance_weights,"list")){
    if ((inherits(distance_weights[[1]],"character") & grepl("\\.tif$",distance_weights[[1]])) |
        inherits(distance_weights[[1]],"PackedSpatRaster")) {
      distance_weights<-lapply(distance_weights,terra::rast)
    } else if (inherits(distance_weights[[1]],"SpatRaster")) {
      distance_weights<-distance_weights
    }
  } else if (inherits(distance_weights,"character")){
    if (grepl("\\.tif$",distance_weights)) {
      distance_weights<-lapply(distance_weights,terra::rast)
    } else if (grepl("\\.rds$",distance_weights)) {
      distance_weights<-readRDS(distance_weights)
    } else {
      if (grepl("\\.zip$",distance_weights)) {
        fls<-unzip(distance_weights,list=T)
        fls<-file.path("/vsizip",distance_weights,fls$Name)
        distance_weights<-lapply(fls,terra::rast)
        names(distance_weights)<-sapply(distance_weights,names)
      } else stop("'distance_weights' must be specified as .rds or .zip file, or a list of file paths or SpatRast objects")
    }
  }

  distance_weights<-lapply(distance_weights,function(x) if (inherits(x,"PackedSpatRaster")) terra::rast(x) else x)
  names(distance_weights)<-lapply(distance_weights,names)

  # Prepare roi data ------------------------------------------------------
  if (is.null(roi)) {
    roi<-distance_weights[[1]]
    roi[!is.na(roi)]<-1
    names(roi)<-"roi"
  } else {
    roi<-process_input(input = roi,
                       input_name="roi",
                       #target = distance_weights[[1]],
                       clip_region=distance_weights[[1]],
                       resample_type="near",
                       working_dir=own_tempdir)
    if (length(roi)>1 & inherits(roi,"SpatRaster")) roi<-roi[[1]]
    if (inherits(roi,"SpatVector") && ncol(roi)>1) roi<-roi[,names(roi)[1]]
    names(roi)<-"roi"
  }

  remove_region<-process_input(input = remove_region, #PS: I'm debating if this is necessary in this function, wouldn't the workflow make more sense to do this external to the function, and feed your desired region into hydroweight()?
                               input_name="remove_region",
                               target = distance_weights[[1]],
                               clip_region=distance_weights[[1]],
                               resample_type="near",
                               working_dir=own_tempdir)

  if (!is.null(remove_region)) {
    roi<-process_input(input = roi,
                       input_name="roi",
                       target = distance_weights[[1]],
                       clip_region=distance_weights[[1]],
                       resample_type="near",
                       working_dir=own_tempdir)
    roi<-roi[[1]]
    names(roi)<-"roi"

    remove_region<-remove_region[[1]]
    names(remove_region)<-"remove_region"

    roi <- terra::mask(x=roi, mask = remove_region, inverse = TRUE)
  }

  roi<-process_input(input = roi, # this converts ROI to polygon
                     input_name="roi",
                     target = terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=terra::crs( distance_weights[[1]])),
                     resample_type="near")

  # Prepare loi -------------------------------------------------------------
  loi<-process_input(input = loi,
                     input_name="loi",
                     variable_names=loi_columns,
                     target = distance_weights[[1]],
                     clip_region=roi,
                     resample_type=loi_resample,
                     working_dir=own_tempdir)

  distance_weights<-lapply(distance_weights,process_input,target=distance_weights[[1]],clip_region=roi,resample_type="bilinear",working_dir=own_tempdir)

  loi_stats <- list(setNames(roi_uid,roi_uid_col) %>% unlist())
  names(loi_stats)<-roi_uid_col
  loi_stats<-list(UID=loi_stats)


  # For raster data with loi_numeric numbers --------------------------------
  if (loi_numeric == TRUE){
    # Non-weighted statistics -------------------------------------------------

    if (any(loi_numeric_stats %in% "mean")){
      loi_mean <- terra::global(loi, fun = "mean", na.rm = T) %>% unlist()
      names(loi_mean)<-paste0(names(loi), "_lumped_mean")
    } else {
      loi_mean<-NULL
    }
    if (any(loi_numeric_stats %in% "sd")){
      loi_sd <- terra::global(loi, fun = "sd", na.rm = T) %>% unlist()
      names(loi_sd)<-paste0(names(loi), "_lumped_sd")
    } else {
      loi_sd<-NULL
    }
    if (any(loi_numeric_stats %in% "median")){
      loi_median <- terra::global(loi, fun = function(x) median(x,na.rm=T)) %>% unlist()
      names(loi_median)<-paste0(names(loi), "_lumped_median")
    }else {
      loi_median<-NULL
    }
    if (any(loi_numeric_stats %in% "min")){
      loi_min <- terra::global(loi, fun = "min", na.rm = T) %>% unlist()
      names(loi_min)<-paste0(names(loi), "_lumped_min")
    }else {
      loi_min<-NULL
    }
    if (any(loi_numeric_stats %in% "max")){
      loi_max <- terra::global(loi, fun = "max", na.rm = T) %>% unlist()
      names(loi_max)<-paste0(names(loi), "_lumped_max")
    }else {
      loi_max<-NULL
    }
    if (any(loi_numeric_stats %in% "sum")){
      loi_sum <- terra::global(loi, fun = "sum", na.rm = T) %>% unlist()
      names(loi_sum)<-paste0(names(loi), "_lumped_sum")
    }else {
      loi_sum<-NULL
    }
    if (any(loi_numeric_stats %in% "cell_count")){
      loi_cell_count <- !is.na(loi)
      loi_cell_count <- terra::global(loi_cell_count, fun="sum", na.rm = T) %>% unlist()
      names(loi_cell_count)<-paste0(names(loi), "_lumped_cell_count")
    }else {
      loi_cell_count<-NULL
    }
    if (any(loi_numeric_stats %in% "NA_cell_count")){
      roi_cell_count <- terra::global(roi, fun="sum", na.rm = T) %>% unlist()
      loi_NA_cell_count <- unlist(roi_cell_count) - unlist(loi_cell_count)
      names(loi_NA_cell_count)<-paste0(names(loi), "_lumped_NA_cell_count")
    }else {
      loi_NA_cell_count<-NULL
    }

    loi_stats <- c(loi_stats,
                   list(lummped=list(mean=loi_mean,
                                     sd=loi_sd,
                                     median=loi_median,
                                     min=loi_min,
                                     max=loi_max,
                                     sum=loi_sum,
                                     cell_count=loi_cell_count,
                                     NA_cell_count=loi_NA_cell_count
                   ))
    )

    # Weighted statistics -----------------------------------------------------
    distance_weights_attributes <- lapply(distance_weights[sapply(distance_weights,names)!="lumped"], function(dw) {

      loi_dist <- loi * dw
      names(loi_dist) <- names(loi)

      loi_sum<-terra::global(loi_dist, fun = "sum", na.rm = TRUE)
      distance_weights_sum<-terra::global(dw, fun = "sum", na.rm = T)

      if (loi_numeric == TRUE) {

        ## Weighted mean
        if (any(loi_numeric_stats %in% c("distwtd_mean","distwtd_sd"))){
          loi_distwtd_mean <- unlist(loi_sum) / unlist(distance_weights_sum)
          names(loi_distwtd_mean)<-paste0(names(loi),"_",names(dw),"_distwtd_mean")
        } else {
          loi_distwtd_mean<-NULL
        }

        if (any(loi_numeric_stats %in% "distwtd_sd")){
          ## Weighted standard deviation
          ## https://stats.stackexchange.com/questions/6534/how-do-i-calculate-a-weighted-standard-deviation-in-excel
          term1 <- terra::global((dw * (loi - unlist(loi_distwtd_mean))^2), fun = "sum", na.rm = TRUE)
          M <- terra::global(dw != 0, fun="sum", na.rm = T)
          term2 <- ((M - 1) / M) * distance_weights_sum

          loi_distwtd_sd <- sqrt(unlist(term1) / unlist(term2))
          names(loi_distwtd_sd)<-paste0(names(loi),"_",names(dw),"_distwtd_sd")
        }else  {
          loi_distwtd_sd<-NULL
        }

        if (!any(loi_numeric_stats %in% c("distwtd_mean"))) loi_distwtd_mean<-NULL
        loi_stats_temp<-list(distwtd_mean=loi_distwtd_mean,
                             distwtd_sd=loi_distwtd_sd)

        if (return_products == TRUE) {
          loi_stats_temp<-c(loi_stats_temp,
                            list(loi_dist_rast=loi_dist))
        }

        return(loi_stats_temp)
      }
    })

    # Add lumped raster to output in return_products=T
    if ("lumped" %in% sapply(distance_weights,names)){
      distance_weights_attributes$lumped<-list(loi_dist_rast=loi)
    }

  }

  # For raster data with loi_numeric categorical numbers -----------------------
  if (loi_numeric == FALSE){

    distance_weights_attributes <- lapply(distance_weights, function(dw) {

      loi_dist <- loi * dw
      names(loi_dist) <- names(loi)

      loi_sum<-terra::global(loi_dist, fun = "sum", na.rm = TRUE)
      distance_weights_sum<-terra::global(dw, fun = "sum", na.rm = T)

      loi_pct_distwtd <- unlist(loi_sum) /  unlist(distance_weights_sum)
      names(loi_pct_distwtd)<-paste0(names(loi),"_",names(dw),"_prop")

      loi_stats_temp<-list(pct_distwtd=loi_pct_distwtd)

      if (return_products == TRUE) {
        loi_stats_temp<-c(loi_stats_temp,
                          list(loi_dist_rast=loi_dist))
      }
      return(loi_stats_temp)
    })
  }



  # Process final table -----------------------------------------------------
  final_out_table<-lapply(distance_weights_attributes,function(x) x[names(x) != "loi_dist_rast"])
  final_out_table<-c(loi_stats,final_out_table)
  final_out_table_nms<-lapply(names(final_out_table),function(x) lapply(names(final_out_table[[x]]),function(y) names(final_out_table[[x]][[y]])))
  final_out_table_nms<-unlist(final_out_table_nms)

  final_out_table<-tibble::enframe(unlist(final_out_table)) %>%
    tidyr::pivot_wider() %>%
    mutate(across(c(everything(),-any_of(roi_uid_col)),as.numeric)) %>%
    mutate(across(ends_with("_prop"),~ifelse(is.na(.),0,.)))

  colnames(final_out_table)<-final_out_table_nms

  final_out_table<-distinct(final_out_table)

  return_list <- list(attribute_table=final_out_table)

  if (return_products == TRUE) {
    final_out_rasts<-lapply(distance_weights_attributes,function(x) x[names(x) == "loi_dist_rast"])
    final_out_rasts<-lapply(final_out_rasts,function(x) lapply(x,terra::wrap))

    return_list<-c(return_list,
                   list(return_products=final_out_rasts))
  }


  # Cleanup -----------------------------------------------------------------
  list_obj<-ls()
  list_obj<-list_obj[!list_obj %in% c("own_tempdir","return_list")]
  rm(list_obj)

  file.remove(list.files(own_tempdir,pattern="_TEMP.tif"))
  fls_to_remove<-unlink(own_tempdir,recursive=T,force =T)
  terra::tmpFiles(current = T,orphan=T,old=F,remove = T)

  return(return_list)
}

