test_that("two people living in close proximity with lots of mosquitoes can be infected by each other.", {
  set.seed(0)
  houses = tibble(x = c(0.4, 0.6), y = 0.5)
  people = tibble(
    ID = 1:2,
    location_ix = list(1,2),
    location_proportions = list(1, 1),
    t_infection = c(0, NA)
  )
  
  # Construct extremely infectious mosquito environment
  mosquito_raster = raster(
    vals = 99999,
    ncols = 1,
    nrows = 1,
    xmn = 0,
    xmx = 1,
    ymn = 0,
    ymx = 1)
  
  sim = Simulation$new(population = people,
                       locations = houses,
                       mosquito = Mosquito$new(bite_rate = 1e5),
                       mosquito_raster = mosquito_raster)
  
  # Advance time
  for (i in 1:30) sim$iterate(1)
  
  # At least one new person has been infected
  expect_equal(sum(sim$linelist$t_infection > 0) > 0, TRUE)
  # At least one person's source is Transmission
  expect_equal("Transmission" %in% sim$linelist$source, TRUE)
})
