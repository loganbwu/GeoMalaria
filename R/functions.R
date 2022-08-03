my_extract = function(...) {
  values = raster::extract(...)
  ifelse(is.null(values), NA_real_, values)
}
