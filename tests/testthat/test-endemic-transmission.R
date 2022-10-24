test_that("under reasonable circumstances transmission occurs indefinitely", {
  set.seed(0)
  n = 100
  rads = seq(0, 2*pi, length.out=n)
  houses = tibble(x = 0.4*cos(rads), y = 0.4*sin(rads))
  people = tibble(
    ID = 1:n,
    location_ix = as.list(1:n),
    location_proportions = as.list(rep(1, n)),
    t_infection = c(0, rep(NA, n-1))
  )
  
  # Construct mosquito environment
  mosquito_raster = raster(
    vals = 1,
    ncols = 1,
    nrows = 1,
    xmn = -1,
    xmx = 1,
    ymn = -1,
    ymx = 1)
  
  # Set a custom immunity profile with complete gradual waning
  Simulation$set("public", "human_immunity", function(population, t) {
    immunity = approx(c(0, 60, 120),
                      c(1, 1, 0),
                      t - population$t_infection,
                      yleft = 0, yright = 1)$y
    replace_na(immunity, 0)
  }, overwrite = TRUE)
  
  sim_base = Simulation$new(population = people,
                            locations = houses,
                            mosquito = Mosquito$new(bite_rate = 1,
                                                    death_rate = 0.25),
                            mosquito_raster = mosquito_raster,
                            p_relapse = 0.5)
  
  # Make copies of the initialised simulation
  sim_1 = sim_base$clone()
  sim_2 = sim_base$clone()
  
  # Simulate both for a long duration but one in weekly time steps, one daily
  n_days = 700
  for (i in 1:n_days) {
    sim_1$iterate(1)
  }
  for (i in 1:(n_days/7)) {
    sim_2$iterate(7)
  }
  
  # Confirm that the daily simulation had transmission in the last few month
  transmission_last_month_1 = with(sim_1$linelist,
                                   length(t_infection[t_infection >= n_days-30 & source=="Transmission"]))
  expect_gt(transmission_last_month_1, 0)
  
  # Confirm that the weekly simulation had transmission in the last few month
  transmission_last_month_2 = with(sim_2$linelist,
                                   length(t_infection[t_infection >= n_days-30 & source=="Transmission"]))
  expect_gt(transmission_last_month_2, 0)
})
