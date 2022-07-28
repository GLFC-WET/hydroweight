#' Process Input GIS files into target format and clip
#'
#' @export
#'
process_input<-function(input=NULL,
                        input_name=NULL,
                        variable_names=NULL,
                        target=NULL,
                        clip_region=NULL,
                        resample_type=c("bilinear","near"),
                        working_dir=NULL,
                        ...){
  require(terra)
  require(sf)

  if (is.null(working_dir)) {
    tdir<-file.path(tempfile())
  } else {
    tdir<-file.path(working_dir,basename(tempfile()))
  }
  if (!dir.exists(tdir)) dir.create(tdir)

  terra::terraOptions(tempdir = tdir, verbose=F)

  if (is.null(input)) return(input)
  output<-NULL

  resample_type<-match.arg(resample_type)

  # Reconcile Types ---------------------------------------------------------
  if (inherits(input,"character")){
    if (grepl("\\.shp$",input)) {
      output<-terra::vect(input)
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
  if (inherits(input,"sf")){
    output<-terra::vect(input)
  }
  if (inherits(input,c("SpatRaster","SpatVector"))){
    output<-input
  }

  if (is.na(terra::crs(output)) | is.null(terra::crs(output))) {
    stop("'output' crs() is NULL or NA. Apply projection before continuing")
  }

  if (is.null(variable_names)) variable_names<-names(output)
  if (any(!variable_names %in% names(output))) stop("some 'variable_names' not in input")

  if (length(variable_names)>0) {
    if (inherits(output,"SpatVector")) output<-output[,variable_names]
    if (inherits(output,"SpatRaster")) output<-terra::subset(output,variable_names)
  }


  if (!is.null(target)){
    target<-process_input(input=target)

    if (is.na(terra::crs(target)) | is.null(terra::crs(target))) {
      stop("'target' crs() is NULL or NA. Apply projection before continuing")
    }

    # Process Clip Region -----------------------------------------------------
    if (!is.null(clip_region)){
      clip_region<-process_input(clip_region,target=terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=terra::crs(target))) #clipping by polygon is slower than by raster
      #clip_region<-process_input(clip_region,target=target)

      if (inherits(output,"SpatVector")) {
        output<-terra::crop(
          x=output,
          y=terra::ext(clip_region)
        )
      }

      if (inherits(output,"SpatRaster")) {
        output<-terra::crop(
          x=output,
          y=terra::ext(clip_region),
          mask=F,
          overwrite=T
        )
      }
    }

    # If target is raster -----------------------------------------------------
    if (inherits(target,"SpatRaster")){
      if (inherits(output,"SpatVector")) {

        if (resample_type=="near"){ # For categorical vector inputs
          output_split<-lapply(setNames(variable_names,variable_names), function(x) {
            out<-output %>%
              st_as_sf() %>%
              select(any_of(x)) %>%
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
            sv<-terra::writeRaster(out,filename=fl,overwrite=T)
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
            out<-terra::rasterize(x=output,
                                  y=target,
                                  field = x,
                                  overwrite=T,
                                  ...
            )
            sv<-terra::writeRaster(out,filename=fl,overwrite=T)
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
        y=clip_region,
        mask=T,
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
    final_temp<-paste0(tempfile(),".tif")
    terra::writeRaster(output,final_temp,overwrite=T)
  }
  if (inherits(output,"SpatVector")){
    final_temp<-paste0(tempfile(),".shp")
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
