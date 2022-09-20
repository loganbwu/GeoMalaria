test_that("repeated relapse behaviour is as expected", {
  set.seed(0)
  # Use Ross Macdonald assumptions
  
  n = 1000
  houses = tibble(x = 0, y = 0)
  people = tibble(
    ID = 1:n,
    location_ix = as.list(rep(1, n)),
    location_proportions = as.list(rep(1, n)),
    t_infection = NA,
    t_relapse = 1
  )
  
  # Construct mosquito environment with NO transmission
  mosquito_raster = raster(
    vals = 0,
    ncols = 1,
    nrows = 1,
    xmn = -1,
    xmx = 1,
    ymn = -1,
    ymx = 1)
  
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 0.1,
                       bite_rate = 0.1,
                       mosquito_death_rate = 0.5,
                       p_relapse = 0.5,
                       mean_recovery = 14)
  
  # Simulate both for the same total time period
  end_t = 25
  diff = pbapply::pbsapply(1:100, function(x) {
    sim_0.5 = sim$clone()
    sim_1 = sim$clone()
    while (sim_0.5$t <= end_t) {
      sim_0.5$iterate(0.75)
    }
    while (sim_1$t <= end_t) {
      sim_1$iterate(1)
    }
    nrow(sim_1$linelist) - nrow(sim_0.5$linelist)
  })
  diff_interval = quantile(diff, c(0.025, 0.975))
  
  # Expect most of the differences in outbreak size are around 0
  expect_true(diff_interval[1] <= 0 & diff_interval[2] >= 0)
})
