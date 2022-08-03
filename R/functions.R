#' Plot spread of mosquitoes over distance
plot_mosquito_migration = function() {
  DX = seq(-1.1 * max_mosquito_flight_range,
           1.1 * max_mosquito_flight_range,
           length.out = 100)
  DT = seq(0, max_mosquito_lifespan, length.out=8)
  DXDT = expand.grid(dx=DX, dt=DT[-1]) %>%
    mutate(density = mosquito_migration(dx, dt))
  
  ggplot(DXDT, aes(x=dx, y=density, color=dt, group=dt)) +
    geom_line() +
    geom_vline(xintercept=c(-1,1)*max_mosquito_flight_range, linetype="dashed")
}


#' Lifecycle stages for P. vivax in humans after inoculation
#'
#' @param t days since inoculation
#' 
#' @return vector 1 if human is infectious, 0 if not
# human_infectivity = function(dt, max_human_infectivity) {
#   t_incubation = runif(length(dt), 12, 17)
#   infectivity = as.numeric(dt > t_incubation & dt < max_human_infectivity)
#   replace_na(infectivity, 0)
# }

my_extract = function(...) {
  values = raster::extract(...)
  ifelse(is.null(values), NA_real_, values)
}
