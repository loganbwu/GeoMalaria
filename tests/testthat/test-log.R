test_that("logging works", {
  set.seed(0)
  n = 100
  rads = seq(0, 2*pi, length.out=n)
  houses = tibble(X = rnorm(n, 0, 15),
                  Y = rnorm(n, 0, 15))
  people = tibble(
    ID = 1:n,
    location_ix = as.list(1:n),
    location_proportions = as.list(rep(1, n)),
    t_infection = c(0, rep(NA, n-1))
  )
  
  # Construct mosquito environment
  mosquito_raster = make_perlin_mosquitoes(list(xmin = floor(min(houses$X)),
                                                xmax = ceiling(max(houses$X)),
                                                ymin = floor(min(houses$Y)),
                                                ymax = ceiling(max(houses$Y))))
  
  # Create two simulations
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 60,
                       bite_rate = 1,
                       mosquito_death_rate = 0.25,
                       p_relapse = 0.1,
                       log_options = c("linelist", "compartment", "EIR"))
  
  # Simulate transmission for a while
  sim$iterate(10)
  
  niter = 5
  dt = 2.5
  days = c(0, 10)
  for (i in 1:niter) {
    sim$iterate(dt)
    days = c(days, sim$t)
  }
  
  # Check linelist data
  expect_identical(sort(names(sim$linelist)), c("ID", "source", "t_infection"))
  expect_gte(nrow(sim$linelist), 1)
  
  # Check compartments data
  expect_identical(sort(names(sim$compartment)), c("infected", "t"))
  expect_identical(sim$compartment$t, days)
  
  # Check animation data
  t_table = table(sim$EIR$t)
  expect_equal(length(unique(t_table)), 1) # All times have the same number of EIR values
  expect_equal(names(t_table), as.character(days)) # Entries for every time step
  expect_type(sim$EIR$ento_inoculation_rate, "double")
  
  # Test a sim with no logging enabled
  sim_nolog = Simulation$new(humans = people,
                             locations = houses,
                             mosquito_raster = mosquito_raster,
                             duration_human_infectivity = 60,
                             bite_rate = 1,
                             mosquito_death_rate = 0.25,
                             p_relapse = 0.1,
                             log_options = NULL)
  sim_nolog$iterate(1)
  expect_identical(sim_nolog$log, list())
  
  
  # Potentially split the actual animation test out
  skip_on_ci() # Don't run on Github Actions
  file <- withr::local_tempfile(
    fileext = ".mkv"
  )
  plot_anim(sim, file)
  expect_true(file.exists(file))
})
