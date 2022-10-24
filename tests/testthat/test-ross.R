test_that("replicate stochastic Ross-Macdonald behaviour", {
  skip("Test not complete")
  set.seed(0)
  # Use Ross Macdonald assumptions
  
  H = 10000
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
  
  # Ross-Macdonald parameters
  rm = list(
    m = M/H,
    a = 0.1,
    b = 1,
    c = 0.1,
    p = 0.75, # P(survive 1 day). g = -log(p)
    r = 0.1, # recovery rate
    # no_recovery = exp(-0.1), # exp(-r), this proportion of people don't recover each day
    v = 0
  )
  R0 = with(rm, {
    (m * a^2 * b * c) / (-log(p) * r) * p^v
  })
  
  # Person is immune if they have a current infection
  Simulation$set("public", "human_immunity", function(population, t) {
    immunity = as.numeric(!is.na(population$t_infection))
    replace_na(immunity, 0)
  }, overwrite = TRUE)
  # mean_infectivity = 1 / rm$r
  
  # Remove mosquito migration
  Simulation$set("public", "mosquito_migration", function(dx, dy, dt) {
    1
  }, overwrite = TRUE)
  
  # Make human infectivity constant
  Simulation$set("public", "human_infectivity", function(population, t) {
    infectivity = ifelse(t - population$t_infection >= 0, 1, 0)
    replace_na(infectivity, 0)
  }, overwrite = TRUE)
  
  # Instant sporogony
  Simulation$set("public", "mosquito_sporogony", function(dt) {
    infectivity = ifelse(dt >= 0, 1, 0)
    replace_na(infectivity, 0)
  }, overwrite = TRUE)
  
  sim = Simulation$new(population = people,
                       locations = houses,
                       mosquito = Mosquito$new(bite_rate = rm$a,
                                               death_rate = 1 - rm$p),
                       mosquito_raster = mosquito_raster,
                       p_relapse = 0,
                       mean_recovery = 14)
  
  
  # Simulate both for a long duration but one in weekly time steps, one daily
  n_days = 10
  pb = progress::progress_bar$new(total = n_days)
  for (i in 1:n_days) {
    pb$tick()
    sim$iterate(0.5)
  }
  plot_epicurve(sim)
  
  Reff = sim$compartment %>%
    mutate(Reff = infected / lag(infected))
  with(Reff, qplot(t, Reff))
  with(sim$history_states, qplot(t, infected))
})
