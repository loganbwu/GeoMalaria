test_that("overall behaviour is not dependent on time step", {
  skip("Test not complete")
  set.seed(0)
  # Use Ross Macdonald assumptions
  
  H = 1000
  M = 10
  X = round(H/100)
  houses = tibble(x = 0, y = 0)
  people = tibble(
    ID = 1:H,
    location_ix = as.list(rep(1, H)),
    location_proportions = as.list(rep(1, H)),
    t_infection = c(rep(0, X), rep(NA, H-X))
  )
  
  # Construct mosquito environment
  mosquito_raster = raster(
    vals = M,
    ncols = 1,
    nrows = 1,
    xmn = -1,
    xmx = 1,
    ymn = -1,
    ymx = 1)
  
  # Person is immune if they have a current infection
  Simulation$set("public", "human_immunity", function(humans, t) {
    immunity = as.numeric(!is.na(humans$t_infection))
    replace_na(immunity, 0)
  }, overwrite = TRUE)
  
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 100,
                       bite_rate = 0.1,
                       mosquito_death_rate = 0.25,
                       p_relapse = 0,
                       mean_recovery = 14)
  
  sim_0.5 = sim$clone()
  sim_1 = sim$clone
  
  # Simulate both for the same total time period
  end_t = 50
  while (sim_0.5$t <= end_t) {
    sim_0.5$iterate(0.5)
  }
  while (sim_1$t <= end_t) {
    sim_0.5$iterate(1)
  }
  
  plot(sim_0.5)
  plot(sim_1)
})
