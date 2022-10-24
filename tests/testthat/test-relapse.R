test_that("one person can relapse if scheduled manually and relapse if scheduled by an infection", {
  set.seed(0)
  houses = tibble(x = 0.5, y = 0.5)
  people = tibble(
    ID = 1,
    location_ix = list(1),
    location_proportions = list(1),
    t_infection = c(0),
    t_relapse = 5 # Schedule a rather early relapse
  )
  
  # Mosquitoes not needed for relapse
  mosquito_raster = raster(
    vals = 0,
    ncols = 1,
    nrows = 1,
    xmn = 0,
    xmx = 1,
    ymn = 0,
    ymx = 1)
  
  sim = Simulation$new(population = people,
                       locations = houses,
                       mosquito = Mosquito$new(bite_rate = 1e5),
                       mosquito_raster = mosquito_raster,
                       p_relapse = 1)
  
  # Advance simulation for period that includes the first relapse
  for (i in 1:10) sim$iterate(1)
  # Advance a lot to guarantee person is due to relapse again
  sim$iterate(200)
  # Advance to evalute relapse
  sim$iterate(1)
  
  # Person relapses after the initial period then has another relapse scheduled
  expect_equal(sum(sim$linelist$source == "Relapse"), 2)
})
