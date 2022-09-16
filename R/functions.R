#' Extract values from Raster objects
#' 
#' See `raster::extract` but returns NA if there are no values to be extracted.
#' 
#' @param ... Arguments to be passed to `raster::extract`
my_extract = function(...) {
  values = raster::extract(...)
  if (is.null(values)) return(NA_real_)
  else return(values)
}


#' Make perlin raster of mosquitoes
#' 
#' Creates a raster to specified dimensions and resolution containing a layer
#' of perlin noise.
#' 
#' @param env_dimensions List of xmn, xmx, ymn, and ymx
#' @param resolution Optional, cells will have edge lengths of this size
make_perlin_mosquitoes = function(env_dimensions, resolution=1) {
  eps = 1e-5
  xmn = floor((env_dimensions$xmin-eps)/resolution) * resolution
  xmx = ceiling((env_dimensions$xmax+eps)/resolution) * resolution
  ymn = floor((env_dimensions$ymin-eps)/resolution) * resolution
  ymx = ceiling((env_dimensions$ymax+eps)/resolution) * resolution
  mosquito_raster = raster(
    vals = 0,
    ncols = round((xmx-xmn)/resolution),
    nrows = round((ymx-ymn)/resolution),
    xmn = xmn,
    xmx = xmx,
    ymn = ymn,
    ymx = ymx)
  noise = noise_perlin(dim(mosquito_raster)[1:2])
  mosquito_raster[] = (noise - min(noise)) * 10000
  
  mosquito_raster
}
