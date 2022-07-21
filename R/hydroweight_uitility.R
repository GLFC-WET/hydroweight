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
                        ...){
  require(terra)
  require(sf)

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
  if (inherits(input,"SpatRaster")){
    output<-input
  }
  if (inherits(input,"SpatVector")){
    output<-input
  }

  if (is.na(terra::crs(output)) | is.null(terra::crs(output))) {
    stop("'output' crs() is NULL or NA. Apply projection before continuing")
  }

  if (is.null(variable_names)) variable_names<-names(output)
  if (any(!variable_names %in% names(output))) stop("some 'variable_names' not in input")

  if (length(variable_names)>0) {
    if (inherits(output,"SpatVector")) output<-output[names(output) %in% variable_names]
    if (inherits(output,"SpatRaster")) output<-terra::subset(output,variable_names)
  }

  if (!is.null(target)){
    target<-process_input(input=target)

    if (is.na(terra::crs(target)) | is.null(terra::crs(target))) {
      stop("'target' crs() is NULL or NA. Apply projection before continuing")
    }

    # If target is raster -----------------------------------------------------
    if (inherits(target,"SpatRaster")){
      if (inherits(output,"SpatVector")) {

        if (resample_type=="near"){ # For categorical vector inputs
          output_split<-lapply(setNames(variable_names,variable_names), function(x) {
            output %>%
              st_as_sf() %>%
              select(any_of(x)) %>%
              split(.[[x]]) %>%
              lapply(vect)
          })

          output_split_nms<-sapply(names(output_split),function(x) paste0(x,"_",names(output_split[[x]])))

          output_split<-unlist(output_split,recursive = F)
          names(output_split)<-output_split_nms

          output <- lapply(output_split, function(x) {
            terra::rasterize(x=x,
                             y=target,
                             field = "",
                             overwrite=T,
                             ...
            )
          })
        }
        if (resample_type=="bilinear"){ # For numeric vector inputs
          output <- lapply(variable_names, function(x) {
            terra::rasterize(x=output,
                             y=target,
                             field = x,
                             overwrite=T,
                             ...
            )
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

    if (inherits(output,"SpatRaster")) {
      if (length(resample_type)>1) stop("'resample_type' must be one of: 'bilinear','near'")
      output <- terra::project(
        x = output,
        y = target,
        method = resample_type,
        overwrite = TRUE,
        ...
      )
    }

    if (inherits(output,"SpatVector")) {
      output <- terra::project(
        x = output,
        y = target,
        ...
      )
    }

  }

  if (!is.null(clip_region) & !is.null(target)){
    clip_region<-process_input(clip_region,target=terra::vect("POLYGON ((0 -5, 10 0, 10 -10, 0 -5))",crs=crs(target)))

    output<-terra::crop(
      x=output,
      y=clip_region,
      mask=T
    )
  }

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
      return(out)
    })

    output <- terra::rast(brick_list)
  }

  return(output)
}
