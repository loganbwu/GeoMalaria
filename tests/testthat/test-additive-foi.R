test_that("Rate of infection is linear with local infectors", {
  set.seed(0)
  
  # Setup for all
  houses = tibble(x = 0, y = 0)
  mosquito_raster = raster(
    vals = 1000,
    ncols = 1,
    nrows = 1,
    xmn = -1,
    xmx = 1,
    ymn = -1,
    ymx = 1)
  max_t = 40 # tune this to get an attack rate where we can discern a difference
  iterations = 100
  
  # Experiment with 1 infector
  n_infected = 1
  n = n_infected + 1
  people = tibble(
    ID = 1:n,
    location_ix = as.list(rep(1, n)),
    location_proportions = as.list(rep(1, n)),
    t_infection = c(rep(0, n_infected), rep(NA, n-n_infected))
  )
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 100,
                       bite_rate = 0.1,
                       mosquito_death_rate = 0.25,
                       p_relapse = 0,
                       mean_recovery = 14)

  rep_results = sapply(seq_len(iterations), function(x) {
    sim_clone = sim$clone()
    sim_clone$iterate(10)
    while (sim_clone$t < max_t) {
      sim_clone$iterate(1)
    }
    # Return whether the nth person was infected
    n %in% sim_clone$linelist$ID
  })
  p_infection = mean(rep_results)
  
  # Experiment with 2 infectors
  n_infected = 2
  n = n_infected + 1
  people = tibble(
    ID = 1:n,
    location_ix = as.list(rep(1, n)),
    location_proportions = as.list(rep(1, n)),
    t_infection = c(rep(0, n_infected), rep(NA, n-n_infected))
  )
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 100,
                       bite_rate = 0.1,
                       mosquito_death_rate = 0.25,
                       p_relapse = 0,
                       mean_recovery = 14)
  
  rep_results = sapply(seq_len(iterations), function(x) {
    sim_clone = sim$clone()
    sim_clone$iterate(10)
    while (sim_clone$t < max_t) {
      sim_clone$iterate(1)
    }
    n %in% sim_clone$linelist$ID
  })
  expectation = 1 - (1 - p_infection) ^ 2
  expectation_sim = rbinom
  actual = mean(rep_results)
  # print(paste0("Expected ", expectation, ", got ", actual))
  
  # Check that the rate of infection with two infectors is approx double the rate of one infector
  # P_2 = 1 - (1 - P_1)^2
  expect_equal(expectation, actual, tolerance = 0.1)
})
