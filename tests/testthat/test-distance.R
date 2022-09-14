test_that("a person nearby can be infected but a person far away cannot be.", {
  set.seed(0)
  
  # Contrived setup where one susceptible person is near and one is very far
  houses = tibble(X = c(0.4, 0.6, 99.5), Y = 0.5)
  people = tibble(
    ID = 1:3,
    location_ix = list(1, 2, 3),
    location_proportions = list(1, 1, 1),
    t_infection = c(0, NA, NA)
  )
  
  # Construct extremely infectious mosquito environment
  mosquito_raster = raster(
    vals = 99999,
    ncols = 1,
    nrows = 1,
    xmn = 0,
    xmx = 100,
    ymn = 0,
    ymx = 1)
  
  sim = Simulation$new(humans = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       duration_human_infectivity = 60,
                       bite_rate = 99999,
                       mosquito_death_rate = 0.25)
  
  # Advance time
  for (i in 1:10) sim$iterate(5)
  
  # Only person #2 has been infected
  with(sim$history_infections,
       expect_equal(
         source[ID == 2] == "Transmission" & all(ID <= 2),
         TRUE)
  )
})