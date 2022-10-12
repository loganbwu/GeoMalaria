test_that("overall behaviour is not dependent on time step size", {
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
  
  ####
  # Setup for first experiment
  n = 1
  n_infected = ceiling(n/100)
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
  
  sim_A = sim$clone()
  sim_B = sim$clone()
  
  # Advance sims in different time steps
  end_t = 50
  while (sim_A$t < end_t) {
    sim_A$iterate(0.5)
  }
  while (sim_B$t < end_t) {
    sim_B$iterate(2)
  }
  # Check that the one reaches identical states
  expect_equal(sim_A$humans, sim_B$humans)
  
  ####
  # Setup for second experiment
  set.seed(0)
  n = 1000
  n_infected = 1 #ceiling(n/100)
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
  sim$iterate(10)
  
  # Advance sims in different time steps
  end_t = 40
  diff = pbapply::pbsapply(1:100, function(x) {
    sim_A = sim$clone()
    sim_B = sim$clone()
    while (sim_A$t < end_t) {
      sim_A$iterate(0.5)
    }
    while (sim_B$t < end_t) {
      sim_B$iterate(2)
    }
    c(A = sim_A$attack_rate, B = sim_B$attack_rate)
  }) %>% t() %>% as_tibble() %>%
    pivot_longer(cols=c("A", "B"), names_to = "Sim", values_to = "attack_rate")
  
  diff_stat = diff %>%
    group_by(Sim) %>%
    summarise(mean = mean(attack_rate),
              var = var(attack_rate),
              sd = sd(attack_rate))
  # Ensure the attack rates after 30 days are approximately equal of time step size
  # Note this is not a t-test. The means will almost surely not be from the same distribution
  expect_true(with(diff_stat,
                   abs(mean[1] - mean[2]) < 1.96 * sqrt(sum(var))))
})
