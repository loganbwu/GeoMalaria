test_that("distribution of relapse delays doesn't depend on which number relapse it is", {
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
  
  p_relapse = 0.4
  sim = Simulation$new(population = people,
                       locations = houses,
                       mosquito_raster = mosquito_raster,
                       p_relapse = p_relapse,
                       mean_recovery = 14)
  
  for (i in 1:1000) {
    sim$iterate(1)
  }
  
  # Check that relapse times make sense
  relapses = sim$linelist %>%
    group_by(ID) %>%
    mutate(ix = row_number() - 1,
           delay = t_infection - lag(t_infection)) %>%
    slice(-1)
  
  # ggplot(relapses, aes(x=t_infection)) +
  #   geom_density() +
  #   facet_wrap(vars(n))

  # ggplot(relapses, aes(x=delay)) +
  #   geom_density() +
  #   facet_wrap(vars(ix))
  
  # Check the distribution of relapse delays
  distributions = relapses %>%
    group_by(ix) %>%
    summarise(mean_delay = mean(delay),
              sd = sd(delay),
              n = n()) %>%
    filter(n > 1)
  means = with(
    distributions,
    list(mean_delay = sum(mean_delay * n) / sum(n),
         mean_sd = sum(sd * n) / sum(n))
  )
  distributions = distributions[distributions$n > 100,]
  differences = list(
    mean_delay = abs(distributions$mean_delay - means$mean_delay),
    mean_sd = abs(distributions$sd - means$mean_sd))
  # All should be roughly similar
  expect_true(all(differences$mean_delay < 5) & all(differences$mean_sd < 5))
  
  
  # ggplot(relapses, aes(x = n)) +
  #   geom_bar()
  
  n_relapses = tibble(n = table(relapses$ix)) %>%
    mutate(frac = n / lag(n)) %>%
    drop_na()
  # Investigate the number of relapses that result in subsequent relapses
  weighted_p = with(n_relapses, sum(n*frac) / sum(n))
  
  # Should be roughly similar
  expect_equal(weighted_p, p_relapse, tolerance=0.2)
})
