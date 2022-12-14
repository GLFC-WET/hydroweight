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
#' digital elevation model/flow accumulation CRS. _TEMP_* files are generated in
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
#' @param hydroweight_dir character. File path to write resulting \code{*.zip} file.
#' @param target_O character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/sample_point.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Target for iEucO or iFLO.
#' @param target_S character (ffull file path with extension, e.g., "C:/Users/Administrator/Desktop/stream_lines.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Target for iEucS or iFLS.
#' @param target_uid character. Unique identifier to precede exported list \code{*.rds} (i.e., "target_uid"_inv_distances.rds)
#' @param OS_combine logical. Should target_O and target_S be merged as targets for iEucS, iFLS, and/or HAiFLS? Use \code{TRUE} or \code{FALSE}. This allows cells surrounding \code{target_O} to flow directly into \code{target_O} rather than be forced through \code{target_S}.
#' @param clip_region character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/lu.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of ESRI Shapefile type or GeoTiFF type only. Region over which distances are calculated.
#' @param dem character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/dem.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster} of GeoTiFF type. Digital elevation model raster.
#' @param flow_accum  character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/flow_accum.tif"), \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}  of GeoTiFF type. Flow accumulation raster (units: # of cells).
#' @param weighting_scheme character. One or more weighting schemes: c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS")
#' @param inv_function function or named list of functions based on \code{weighting_scheme} names. Inverse function used in \code{terra::app()} to convert distances to inverse distances. Default: \code{(X * 0.001 + 1)^-1} assumes projection is in distance units of m and converts to distance units of km.
#' @param return_products logical. If \code{TRUE}, a list containing the file path to write resulting \code{*.zip} file, and \code{distance_weights} raster. If \code{FALSE}, file path only.
#' @param save_output logical. Should output rasters be saved to a zip file?
#' @param wrap_return_products logical. Should return_products be \code{terra::wrap()} (necessary for running function in parallel)?
#' @param clean_tempfiles logical. Should temporary files be removed?
#'
#' @return Named list of \code{PackedSpatRaster} distance-weighted rasters and location of accompanying \code{*.zip} in \code{hydroweight_dir}
#' @export
#'
hydroweight <- function(hydroweight_dir,
                        target_O = NULL,
                        target_S = NULL,
                        target_uid = NULL,
                        OS_combine = NULL,
                        clip_region = NULL,
                        dem,
                        flow_accum = NULL,
                        weighting_scheme = c("lumped", "iEucO", "iEucS", "iFLO", "iFLS", "HAiFLO", "HAiFLS"),
                        inv_function = function(x) {
                          (x * 0.001 + 1)^-1
                        },
                        return_products = T,
                        wrap_return_products=T,
                        save_output=T,
                        clean_tempfiles=T,
                        temp_dir=NULL) {

  if (!save_output & !return_products) stop("'return_products' and 'save_output' cannot both be FALSE")

  weighting_scheme<-match.arg(weighting_scheme,several.ok = T)

  require(terra)

  if (!is.null(temp_dir)){
    own_tempdir<-tempfile("",temp_dir)
  } else {
    own_tempdir<-gsub("file","",tempfile())
  }

  if (!dir.exists(own_tempdir)) dir.create(own_tempdir)

  ## Set up raster temp file for out-of-memory raster files
  terra::terraOptions(tempdir = own_tempdir, verbose=F)

  ## Set whitebox verbose_mode to FALSE
  whitebox::wbt_options(verbose = FALSE)


  if (is.null(target_uid)) stop("target_uid must be specified")

  # Process DEM Input ----------------------------------------------------------
  dem<-process_input(input = dem,input_name="dem",working_dir=own_tempdir)
  dem<-dem[[1]]

  dem_crs<-terra::crs(dem)

  # Process clip_region Input --------------------------------------------------
  clip_region<-process_input(input = clip_region,
                             #target=dem,
                             #clip_region=dem,
                             target = terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=dem_crs),
                             clip_region = terra::as.polygons(terra::ext(dem),crs=terra::crs(dem)),
                             input_name="clip_region",
                             working_dir=own_tempdir)
  clip_region<-clip_region[,1]

  dem<-process_input(input = dem,
                     target = dem,
                     clip_region=clip_region,
                     resample_type="bilinear",
                     input_name="dem",
                     working_dir=own_tempdir)

  # Process flow_accum Input ---------------------------------------------------
  flow_accum<-process_input(input = flow_accum,
                            target = dem,
                            clip_region=clip_region,
                            resample_type = "bilinear",
                            input_name="flow_accum",
                            working_dir=own_tempdir)
  flow_accum<-flow_accum[[1]]

  # Process target_O Input -----------------------------------------------------
  target_O<-process_input(input = target_O,
                          target = dem,
                          clip_region=clip_region,
                          input_name="target_O",
                          working_dir=own_tempdir)
  target_O<-target_O[[1]]
  target_O[!is.na(target_O)]<-1

  # Process target_S Input -----------------------------------------------------
  target_S<-process_input(input = target_S,
                          target = dem,
                          clip_region=clip_region,
                          input_name="target_S",
                          working_dir=own_tempdir)
  target_S<-target_S[[1]]
  target_S[!is.na(target_S)]<-1

  # Setup OS combined -------------------------------------------------------
  if (is.null(OS_combine)) {
    OS_combine <- FALSE
  }

  if (OS_combine == TRUE & !is.null(target_O) & !is.null(target_S)) {
    (target_S_OS <- target_S)
    (target_O_OS <- target_O)

    target_S_OS[target_S_OS > 0] <- 1
    target_O_OS[target_O_OS > 0] <- 1

    OS_combine_r <- terra::mosaic(x=terra::sprc(list(target_O_OS, target_S_OS)), fun="sum")
    OS_combine_r[OS_combine_r > 0] <- 1
    OS_combine_r[OS_combine_r == 0] <- NA
  } else {
    OS_combine_r<-NULL
  }

  rast_list<-list(
    dem,
    flow_accum,
    clip_region,
    target_O,
    target_S,
    OS_combine_r
  )

  rast_list_nms<-c(
    paste0(target_uid,"_TEMP_dem_clip.tif"),
    paste0(target_uid,"_TEMP_flow_accum_clip.tif"),
    paste0(target_uid,"_TEMP_clip_region.tif"),
    paste0(target_uid,"_TEMP_target_O_clip.tif"),
    paste0(target_uid,"_TEMP_target_S_clip.tif"),
    paste0(target_uid,"_TEMP_OS_combine.tif")
  )
  names(rast_list)<-rast_list_nms

  rast_list_nms<-lapply(rast_list_nms,function(x){
    if (inherits(rast_list[[x]],"SpatVector")) {
      out<-gsub("\\.tif","\\.shp",x)
    } else {
      out<-x
    }
    return(out)
  })
  names(rast_list)<-rast_list_nms

  rast_list<-rast_list[!sapply(rast_list,is.null)]
  rast_list<-rast_list[sapply(rast_list,length)>0]

  # Write temporary rasters -------------------------------------------------
  for (i in names(rast_list)){
    if (!is.null(rast_list[[i]])){
      if (inherits(rast_list[[i]],"SpatVector")) {
        writeVector(rast_list[[i]],file.path(own_tempdir,i),overwrite=T)
      } else {
        writeRaster(rast_list[[i]],file.path(own_tempdir,i),overwrite=T,gdal="COMPRESS=NONE")
      }
    }
  }


  lumped_inv <- dem
  lumped_inv[!is.na(lumped_inv)]<-1
  names(lumped_inv)<-"lumped"
  terra::writeRaster(lumped_inv,file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),overwrite=T,gdal="COMPRESS=NONE")
  terra::writeRaster(lumped_inv,file.path(own_tempdir, paste0(target_uid,"_TEMP_lumped.tif")),overwrite=T,gdal="COMPRESS=NONE")

  ## iEucO, Euclidean distance to target_O ----
  if ("iEucO" %in% weighting_scheme) {

    whitebox::wbt_cost_distance(
      source = file.path(own_tempdir, paste0(target_uid,"_TEMP_target_O_clip.tif")),
      cost = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
      out_accum = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")),
      out_backlink = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_backlink.tif")),
      verbose_mode = FALSE
    )

    iEucO <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")))
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
    iEucO_inv<-terra::mask(iEucO_inv,lumped_inv)
    names(iEucO_inv)<-"iEucO"
    if (save_output) terra::writeRaster(iEucO_inv,
                       file.path(own_tempdir, paste0(target_uid,"_TEMP_iEucO.tif")),
                       overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                       NAflag = -9999
    )
  }

  ## iEucS, Euclidean distance to target_S ----
  if ("iEucS" %in% weighting_scheme) {
    if (OS_combine == FALSE) {
      whitebox::wbt_cost_distance(
        source = file.path(own_tempdir, paste0(target_uid,"_TEMP_target_S_clip.tif")),
        cost = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
        out_accum = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")),
        out_backlink = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_backlink.tif")),
        verbose_mode = FALSE
      )

      iEucS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")))
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
      iEucS_inv<-terra::mask(iEucS_inv,lumped_inv)
      names(iEucS_inv)<-"iEucS"
      if (save_output) terra::writeRaster(iEucS_inv,
                         file.path(own_tempdir, paste0(target_uid,"_TEMP_iEucS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if (OS_combine == TRUE) {
      whitebox::wbt_cost_distance(
        source = file.path(own_tempdir, paste0(target_uid,"_TEMP_OS_combine.tif")),
        cost = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip_cost.tif")),
        out_accum = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")),
        out_backlink = file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_backlink.tif")),
        verbose_mode = FALSE
      )

      iEucS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_cost_distance.tif")))
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
      iEucS_inv<-terra::mask(iEucS_inv,lumped_inv)
      names(iEucS_inv)<-"iEucS"
      if (save_output) terra::writeRaster(iEucS_inv,
                         file.path(own_tempdir, paste0(target_uid,"_TEMP_iEucS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }
  }

  ## iFLO, flow line distance to target_O ----
  if ("iFLO" %in% weighting_scheme | "HAiFLO" %in% weighting_scheme) {
    whitebox::wbt_downslope_distance_to_stream(
      dem = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip.tif")),
      streams = file.path(own_tempdir, paste0(target_uid,"_TEMP_target_O_clip.tif")),
      output = file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLO.tif")),
      verbose_mode = FALSE
    )

    iFLO <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLO.tif")))
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
    names(iFLO_inv)<-"iFLO"
    if (save_output) terra::writeRaster(iFLO_inv,
                       file.path(own_tempdir, paste0(target_uid,"_TEMP_iFLO.tif")),
                       overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                       NAflag = -9999
    )
  }

  ## iFLS, flow line distance to target_S ----
  if ("iFLS" %in% weighting_scheme | "HAiFLS" %in% weighting_scheme) {
    if (OS_combine == FALSE) {
      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip.tif")),
        streams = file.path(own_tempdir, paste0(target_uid,"_TEMP_target_S_clip.tif")),
        output = file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")),
        verbose_mode = FALSE
      )

      iFLS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")))
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
      names(iFLS_inv)<-"iFLS"
      if (save_output) terra::writeRaster(iFLS_inv,
                         file.path(own_tempdir, paste0(target_uid,"_TEMP_iFLS.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if (OS_combine == TRUE) {
      whitebox::wbt_downslope_distance_to_stream(
        dem = file.path(own_tempdir, paste0(target_uid,"_TEMP_dem_clip.tif")),
        streams = file.path(own_tempdir, paste0(target_uid,"_TEMP_OS_combine.tif")),
        output = file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")),
        verbose_mode = FALSE
      )

      iFLS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")))
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
      names(iFLS_inv)<-"iFLS"
      if (save_output) terra::writeRaster(iFLS_inv,
                         file.path(own_tempdir, paste0(target_uid,"_TEMP_iFLS.tif")),
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

    if ("HAiFLO" %in% weighting_scheme) {

      HAiFLO <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLO.tif")))
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

      HAiFLO_inv <- HAiFLO_inv * flow_accum

      terra::crs(HAiFLO_inv) <- dem_crs
      names(HAiFLO_inv)<-"HAiFLO"
      if (save_output) terra::writeRaster(HAiFLO_inv,
                         file.path(own_tempdir, paste0(target_uid,"_TEMP_HAiFLO.tif")),
                         overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                         NAflag = -9999
      )
    }

    if ("HAiFLS" %in% weighting_scheme) {
      if (OS_combine == TRUE) {

        HAiFLS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")))
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

        HAiFLS_inv <- HAiFLS_inv * flow_accum
        HAiFLS_inv <- terra::mask(HAiFLS_inv, OS_combine_r, maskvalues = 1)
        names(HAiFLS_inv)<-"HAiFLS"

        if (save_output) terra::writeRaster(HAiFLS_inv,
                           file.path(own_tempdir, paste0(target_uid,"_TEMP_HAiFLS.tif")),
                           overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                           NAflag = -9999
        )
      }

      if (OS_combine == FALSE) {
        HAiFLS <- terra::rast(file.path(own_tempdir, paste0(target_uid,"_TEMP_flowdist-iFLS.tif")))
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

        HAiFLS_inv <- HAiFLS_inv * flow_accum
        HAiFLS_inv <- terra::mask(HAiFLS_inv, target_S, maskvalues = 1)
        names(HAiFLS_inv)<-"HAiFLS"

        if (save_output) terra::writeRaster(HAiFLS_inv,
                           file.path(own_tempdir, paste0(target_uid,"_TEMP_HAiFLS.tif")),
                           overwrite = TRUE, gdal = c("COMPRESS=NONE"),
                           NAflag = -9999
        )
      }
    }
  }

  ## PREPARE OUTPUT ----

  weighting_scheme_inv <- paste0(weighting_scheme, "_inv")

  dist_list <- lapply(weighting_scheme_inv, function(x) {
    out<-get(x)
    names(out)<-gsub("_inv","",x)
    return(out)
  })
  names(dist_list) <- weighting_scheme

  if (save_output) {
    dist_list_out<-lapply(dist_list,function(x) {
      file.rename(
        file.path(own_tempdir,paste0(paste0(target_uid,"_TEMP_",names(x),".tif"))),
        file.path(own_tempdir,paste0(paste0(names(x),".tif")))
      )
      fl<-file.path(own_tempdir,paste0(paste0(names(x),".tif")))
      #terra::writeRaster(x,fl,overwrite=T,gdal="COMPRESS=NONE")
      return(fl)
    })

    out_file<-file.path(hydroweight_dir,paste0(target_uid, "_inv_distances.zip"))
    if (file.exists(out_file)) {
      out_file<-file.path(hydroweight_dir,paste0(target_uid,"_",basename(tempfile()), "_inv_distances.zip"))
      warning(paste0("Target .zip file already exists. Saving as: ",out_file))
    }
    zip(out_file,
        unlist(dist_list_out),
        flags = '-r9Xjq'
    )
  } else {
    out_file<-NULL
    dist_list_out<-NULL
  }


  if (return_products){
    if (wrap_return_products) {
      dist_list<-lapply(dist_list,terra::wrap) # will need to use wrap() here for terra: https://github.com/rspatial/terra/issues/50 -- this is slow for large rasters
    } else {
      dist_list<-dist_list
    }
  } else {
    dist_list<-out_file
  }

  ## CLEAN UP ----
  if (clean_tempfiles){
    temp_rasters <- list.files(path = own_tempdir, pattern = paste0(target_uid,"_hydroweight_"),
                               full.names = TRUE)
    file.remove(temp_rasters)

    temp_rasters <- list.files(path = own_tempdir, pattern = paste0(target_uid,"_TEMP_"),
                               full.names = TRUE)
    r1<-file.remove(temp_rasters)

    if (!is.null(dist_list_out)){
      r2<-file.remove(unlist(dist_list_out))
    }

    list_obj<-ls()
    list_obj<-list_obj[!list_obj %in% c("own_tempdir","dist_list")]
    rm(list=list_obj)

    fls_to_remove<-unlink(own_tempdir,recursive=T,force =T)

    suppressWarnings(terra::tmpFiles(current = T,orphan=T,old=T,remove = T))

  }

  return(dist_list)
}
