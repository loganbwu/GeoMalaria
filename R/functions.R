my_extract = function(...) {
  values = raster::extract(...)
  if (is.null(values)) return(NA_real_)
  else return(values)
}


# Mosquitoes
make_perlin_mosquitoes = function(env_dimensions, resolution=1) {
  mosquito_raster = raster(
    vals = 0,
    ncols = abs((env_dimensions$xmax - env_dimensions$xmin) / resolution),
    nrows = abs((env_dimensions$ymax - env_dimensions$ymin) / resolution),
    xmn = env_dimensions$xmin,
    xmx = env_dimensions$xmax,
    ymn = env_dimensions$ymin,
    ymx = env_dimensions$ymax,
    crs = "NULL +units=km")
  noise = noise_perlin(dim(mosquito_raster)[1:2])
  mosquito_raster[] = (noise - min(noise)) * 10000
  names(mosquito_raster) = "count"
  
  mosquito_raster
}

