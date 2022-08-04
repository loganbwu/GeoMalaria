my_extract = function(...) {
  values = raster::extract(...)
  if (is.null(values)) return(NA_real_)
  else return(values)
}


# Mosquitoes
make_perlin_mosquitoes = function(env_dimensions) {
  mosquito_raster = raster(vals=0, nrows=100, ncols=100,
                              xmn=0, xmx=env_dimensions[1],
                              ymn=0, ymx=env_dimensions[2],
                              crs=NA)#"+units=km")
  noise = noise_perlin(dim(mosquito_raster)[1:2])
  mosquito_raster[] = (noise - min(noise)) * 10000
  names(mosquito_raster) = "count"
  
  mosquito_raster
}

