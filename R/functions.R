#' Lifecycle stages for P. vivax in mosquitos after a blood meal
#'
#' @param t days since blood meal
#' 
#' @return boolean vector TRUE if mosquito is infectious
vivax_mosquito_lifecycle = function(t) {
  t_sporogony = runif(length(t), 8, 10)
  infectious = t > t_sporogony
  ifelse(is.na(infectious), F, infectious)
}

#' Proportion of mosquitoes that survive over time
mosquito_survival = function(dt) {
  rate = 0.5
  survival = dexp(dt, rate=rate) / rate
}

#' Infectivity of the vivax parasite in mosquitoes over time
mosquito_sporogony = function(dt) {
  t_sporogony = runif(length(dt), 8, 10)
  as.numeric(dt > t_sporogony)
}

#' Overall infectivity of a mosquito after an infected blood meal
mosquito_infectivity = function(dt) {
  mosquito_survival(dt) * mosquito_sporogony(dt)
}


# 1-D variance of mosquito travel in (km/day)^2
mosquito_travel = 0.5^2
mosquito_migration = function(dx, dt) {
  dnorm(dx, sd=sqrt(mosquito_travel*dt))
}

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
#' @return boolean vector TRUE if human is infectious
human_infectivity = function(dt) {
  t_incubation = runif(length(dt), 12, 17)
  infectious = dt > t_incubation & dt < max_human_infectivity
  ifelse(is.na(infectious), F, infectious)
}

my_extract = function(...) {
  values = raster::extract(...)
  ifelse(is.null(values), NA_real_, values)
}
