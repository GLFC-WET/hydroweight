#' Process Input GIS files into target format and clip
#'
#' \code{hydroweight::process_input()} processes any input GIS data files into type \code{SpatVector} or \code{SpatRaster} depending on `target`.
#'
#' `input` data is coerced to type \code{SpatVector} or \code{SpatRaster} depending on class of `target`.
#' If `resample_type` == \code{bilinear}, and `target` is of class \code{SpatRaster}, `input` is
#' coerced to class \code{SpatRaster}, reasterizing and/or reprojecting as necessary. If `resample_type` == \code{near},
#' `target` is of class \code{SpatRaster}, and `input` is of class \code{SpatRaster}, unique values of `input` are coerced
#' into separate layers of a \code{SpatRaster}, with 1 indicating the value is present; if input` is of class \code{SpatVector},
#' unique values of `variable_names` are each coerced into separate layers in a \code{SpatRaster}, with 1 indicating the value is present.
#'
#'
#' @param input character. (full file path with extension, e.g., "C:/Users/Administrator/Desktop/input.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Input GIS data.
#' @param input_name \code{NULL} or character. A name that is used only in error or warning messages.
#' @param variable_names \code{NULL} or character. The names if input that should be included in the coerced output. If \code{NULL} all variables in `input` are used.
#' @param target \code{NULL} or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/target.tif"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Target for intput to be coerced into. If \code{NULL} all variables in `input` are returned as \code{SpatVector} or \code{SpatRaster}.
#' @param clip_region \code{NULL} or character (full file path with extension, e.g., "C:/Users/Administrator/Desktop/clip_region.shp"), \code{sf}, \code{SpatVector}, \code{PackedSpatVector}, \code{RasterLayer}, \code{SpatRaster}, or \code{PackedSpatRaster}. Input with be clipped and masked to this region. Internally converted to  \code{SpatVector}. If \code{NULL}, full extent of `input` is returned.
#' @param resample_type character. If \code{bilinear} (default), the `variable_names` being coerced are numeric, if \code{'near'}, the `variable_names` being coerced are categorical, for the purposes of resampling/projecting if `target` is \code{SpatRaster}.
#' @param working_dir \code{NULL} or character. Folder path for temporary file storage. If \code{NULL}, assigned internally.
#' @param ... other variables passed to various \code{terra} functions.
#' @return an object of class \code{SpatVector} or \code{SpatRaster} depending on `target`
#' @export

process_input<-function(input=NULL,
                        input_name=NULL,
                        variable_names=NULL,
                        target=NULL,
                        clip_region=NULL,
                        resample_type=c("bilinear","near"),
                        working_dir=NULL,
                        ...){
  if (is.null(input)) return(input)
  output<-NULL

  require(terra)
  require(sf)
  require(dplyr)

  if (is.null(working_dir)) {
    tdir<-file.path(gsub("file","",tempfile()))
    working_dir<-tdir
  } else {
    tdir<-file.path(working_dir,basename(gsub("file","",tempfile())))
  }
  if (!dir.exists(tdir)) dir.create(tdir)

  terra::terraOptions(tempdir = tdir, verbose=F)

  resample_type<-match.arg(resample_type)

  # Reconcile Types ---------------------------------------------------------
  if (inherits(input,"character")){
    if (grepl("\\.shp$",input)) {
      output<-terra::vect(input) # terra seems to have trouble reading crs sometimes (sf::read_sf())
    } else if (grepl("\\.tif$|\\.tiff",input)) {
      output<-terra::rast(input)
    } else  {
      stop(paste0(input_name," must be spcified as vector or raster layer, or a character string ending in .shp or .tif"))
    }
  }
  if (inherits(input,"RasterLayer")){
    output<-terra::rast(input)
  }
  if (inherits(input,"PackedSpatRaster")){
    output<-terra::rast(input)
  }
  if (inherits(input,"PackedSpatVector")){
    output<-terra::vect(input)
  }
  if (inherits(input,c("sf","sfc"))){
    output<-terra::vect(input)
  }
  if (inherits(input,c("SpatRaster","SpatVector"))){
    output<-input
  }

  if (is.na(terra::crs(output)) | is.null(terra::crs(output))) {
    stop("'output' crs() is NULL or NA. Apply projection before continuing")
  }

  orig_output<-output

  if (is.null(variable_names)) variable_names<-names(output)
  if (any(!variable_names %in% names(output))) stop("some 'variable_names' not in input")

  if (length(variable_names)>0) {
    if (inherits(output,"SpatVector")) output<-output[,variable_names]
    if (inherits(output,"SpatRaster")) output<-terra::subset(output,variable_names)
  }


  if (!is.null(target)){
    target<-process_input(input=target,working_dir=tdir)
    clip_region<-process_input(clip_region,target=terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=terra::crs(target)),working_dir=tdir)

    if (is.na(terra::crs(target)) | is.null(terra::crs(target))) {
      stop("'target' crs() is NULL or NA. Apply projection before continuing")
    }

    # If target is raster -----------------------------------------------------
    if (inherits(target,"SpatRaster")){
      if (inherits(output,"SpatVector")) {

        output<-terra::project(output,target)

        if (resample_type=="near"){ # For categorical vector inputs
          output_split<-lapply(setNames(variable_names,variable_names), function(x) {
            out<-output %>%
              sf::st_as_sf() %>%
              dplyr::select(tidyselect::any_of(x)) %>%
              split(.[[x]]) %>%
              lapply(terra::vect)

            fl<-lapply(out, function(x) file.path(tdir,paste0(basename(tempfile()),".shp")))
            sv<-lapply(names(out),function(x) terra::writeVector(out[[x]],filename=fl[[x]],overwrite=T))

            out<-lapply(fl,terra::vect)

            out<-lapply(out,function(y) {
              names(y)<-paste0(x,"_",unlist(y[[1]])[[1]])
              return(y)
            })

            return(out)
          })

          output_split_nms<-sapply(names(output_split),function(x) paste0(x,"_",names(output_split[[x]])))

          output_split<-unlist(output_split,recursive = F)
          names(output_split)<-output_split_nms

          output <- lapply(output_split, function(x) {
            fl<-file.path(tdir,paste0(basename(tempfile()),".tif"))
            out<-terra::rasterize(x=x,
                                  y=target,
                                  field = "",
                                  overwrite=T,
                                  ...
            )
            sv<-terra::writeRaster(out,filename=fl,overwrite=T,gdal="COMPRESS=NONE")
            return(terra::rast(fl))
          })

          output<-lapply(setNames(names(output),names(output)),function(x){
            out<-output[[x]]
            names(out)<-x
            return(out)
          })
        }
        if (resample_type=="bilinear"){ # For numeric vector inputs
          output <- lapply(variable_names, function(x) {
            fl<-file.path(tdir,paste0(basename(tempfile()),".tif"))
            out<-terra::rasterize(x=output, #PS: Development version of terra gives a warning here, may cause problems in the future
                                  y=target,
                                  field = x,
                                  overwrite=T,
                                  ...
            )
            sv<-terra::writeRaster(out,filename=fl,overwrite=T,gdal="COMPRESS=NONE")
            return(terra::rast(fl))
          })
        }

        output <- terra::rast(output)
      }
    }

    # If target is vector -----------------------------------------------------
    if (inherits(target,"SpatVector")){
      if (inherits(output,"SpatRaster")) {
        if (terra::is.polygons(target)) output<-terra::as.polygons(output,...)
        if (terra::is.points(target)) output<-terra::as.points(output,...)
        if (terra::is.lines(target)) output<-terra::as.lines(output,...)
      }
    }


    # Reprojecting to mathc target --------------------------------------------
    if (inherits(output,"SpatRaster")) {
      need_reproj<-terra::compareGeom(
        x = output,
        y = target,
        lyrs=F,
        crs=T,
        ext=T,
        rowcol=T,
        res=T,
        warncrs=F,
        stopOnError=F,
        messages=F
      )
      need_reproj<-!need_reproj # take inverse

      if (need_reproj) {
        if (length(resample_type)>1) stop("'resample_type' must be one of: 'bilinear','near'")

        output <- terra::project(
          x = output,
          y = target,
          method = resample_type,
          overwrite = TRUE,
          #filename = file.path(tdir,paste0(basename(tempfile()),".tif")),
          ...
        )
      }

    }

    if (inherits(output,"SpatVector")) {
      output <- terra::project(
        x = output,
        y = target,
        ...
      )
    }

  }

  # Clip and mask to clip_region --------------------------------------------

  if (!is.null(clip_region) & !is.null(target)){
    if (inherits(output,"SpatVector")) {
      output<-terra::crop(
        x=output,
        y=clip_region
      )
    }
    if (inherits(output,"SpatRaster")) {
      output<-terra::crop(
        x=output,
        snap="in",
        y=clip_region,
        mask=T,
        overwrite=T
      )

      output<-terra::mask(
        x=output,
        mask=clip_region,
        overwrite=T
      )
    }
  }


  # Separate output raster into separate layers if multiple numerical values present ------
  if (inherits(output,"SpatRaster")) {
    if (resample_type=="near") {
      if (terra::nlyr(output)>1) {
        brick_list<-terra::split(output,names(output))
        names(brick_list)<-names(output)
      } else {
        brick_list<-output
      }

      brick_list<-lapply(names(brick_list),function(y){
        uv <- unlist(terra::unique(brick_list[[y]]))
        uv <- uv[!sapply(uv,is.nan)]
        uv <- uv[!sapply(uv,is.na)]
        uv <- uv[!sapply(uv,is.infinite)]

        if (length(uv)>1) {
          out<-lapply(uv,function(x) {
            repl<-terra::values(brick_list[[y]])
            repl[repl!=x]<-NA
            repl[!is.na(repl)]<-1
            out<-terra::setValues(brick_list[[y]],repl)
            if (length(uv)>1) names(out)<-paste0(y,"_",x)
            if (length(uv)==1) names(out)<-y
            return(out)
          })

          nms<-sapply(out,names)
          out <- terra::rast(out)
          names(out) <- nms
        } else {
          out<-brick_list[[y]]
        }
        return(out)
      })

      output <- terra::rast(brick_list)
    }
  }

  if (inherits(output,"SpatRaster")){

    if (!inherits(orig_output,"SpatVector")){
      need_save<-terra::compareGeom(
        x = output,
        y = orig_output,
        lyrs=T,
        crs=T,
        ext=T,
        rowcol=T,
        res=T,
        warncrs=F,
        stopOnError=F,
        messages=F
      )
      need_save<-!need_save # take inverse
    } else {
      need_save<-TRUE
    }

    if (need_save) {
      final_temp<-file.path(tdir,paste0(basename(tempfile()),".tif"))
      terra::writeRaster(output,final_temp,overwrite=T,gdal="COMPRESS=NONE")
    }

  }
  if (inherits(output,"SpatVector")){
    final_temp<-file.path(tdir,paste0(basename(tempfile()),".shp"))
    terra::writeVector(output,final_temp,overwrite=T)
  }

  # fls_to_remove<-list.files(tdir,recursive=T,full.names=T)
  # fls_to_remove<-fls_to_remove[!grepl(gsub("\\.tif$|\\.shp$","",final_temp),fls_to_remove)]

  # list_obj<-ls()
  # list_obj<-list_obj[!list_obj %in% c("tdir","output")]
  #
  # rm(list=list_obj)
  # fls_to_remove<-unlink(tdir,recursive=T,force =T)
  # fls_to_remove<-file.remove(fls_to_remove)

  # terra::tmpFiles(current = F,orphan=T,old=F,remove = T)

  return(output)
}
