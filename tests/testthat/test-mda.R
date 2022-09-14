test_that("iff relapse is disabled <=> MDA will stop transmission", {
  set.seed(0)
  n = 100
  rads = seq(0, 2*pi, length.out=n)
  houses = tibble(X = 0.4*cos(rads), Y = 0.4*sin(rads))
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
  
  # Set a custom immunity profile
  complete_waning = function(dt) {
    immunity = approx(c(0, 10),
                      c(1, 1),
                      dt,
                      yleft = 0, yright = 0)$y
    replace_na(immunity, 0)
  }
  Simulation$set("public", "human_immunity", complete_waning, overwrite = TRUE)
  
  # Create two simulations
  sim_norelapse = Simulation$new(humans = people,
                                 locations = houses,
                                 mosquito_raster = mosquito_raster,
                                 duration_human_infectivity = 60,
                                 bite_rate = 1,
                                 mosquito_death_rate = 0.25,
                                 p_relapse = 0)
  
  sim_relapse = Simulation$new(humans = people,
                               locations = houses,
                               mosquito_raster = mosquito_raster,
                               duration_human_infectivity = 60,
                               bite_rate = 1,
                               mosquito_death_rate = 0.25,
                               p_relapse = 0.5)
  
  # Simulate transmission for a while
  for (i in 1:100) {
    sim_norelapse$iterate(1)
    sim_relapse$iterate(1)
  }
  
  # Interrupt and clear all current humand and mosquito infections
  sim_norelapse$humans = sim_norelapse$humans %>%
    mutate(t_infection = NA)
  sim_norelapse$mosquito_infections = sim_norelapse$mosquito_infections %>%
    slice(0)
  
  sim_relapse$humans = sim_relapse$humans %>%
    mutate(t_infection = NA)
  sim_relapse$mosquito_infections = sim_relapse$mosquito_infections %>%
    slice(0)
  
  # Simulate transmission for a while
  for (i in 1:100) {
    sim_norelapse$iterate(1)
    sim_relapse$iterate(1)
  }
  
  # Confirm no infections occur after MDA
  expect_equal(max(sim_norelapse$history_infections$t_infection) <= 100, TRUE)
  # Relapse brings transmission back after MDA
  expect_equal(max(sim_relapse$history_infections$t_infection) > 100, TRUE)
})
