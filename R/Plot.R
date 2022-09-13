print.Simulation = function(sim) {
  print(
    list(
      "t" = sim$t,
      "History" = sim$history_infections
    )
  )
}

plot_init = function(sim) {
  mosquito_data = as.data.frame(sim$mosquito_raster, xy=T) %>%
    rename(X = x, Y = y)
  human_data = sim$humans_expand %>%
    mutate(State = case_when(t_infection == 0 ~ "Importation",
                             TRUE ~ "Susceptible")) %>%
    arrange(t_infection)
  ggplot(mapping = aes(x=X, y=Y)) +
    geom_raster(data = mosquito_data, aes(fill=layer)) +
    geom_point(data = human_data, aes(color=State, size=location_proportions), alpha=0.6) +
    coord_equal() +
    scale_size_area(max_size = 3) +
    scale_fill_viridis_c(option = "magma") +
    scale_color_manual(values = c("Importation"="tomato",
                                  "Susceptible"="grey")) +
    labs(title = "Initial state",
         color = NULL,
         fill = "Mos/km^2",
         size = "Proportion\ntime spent",
         x = NULL,
         y = NULL) +
    theme_minimal() +
    theme(legend.key.height = unit(10, "pt"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())
}

plot.Simulation = function(sim) {
  p_state = plot_state(sim)
  p_epicurve = plot_epicurve(sim)
  p_state / p_epicurve + plot_layout(heights=c(5,1), guides="collect")
}

plot_state = function(sim) {
  humans = sim$humans_collapse %>%
    mutate(Infection = case_when(t_infection == sim$t ~ "New",
                                 ID %in% sim$history_infections$ID ~ "Historical",
                                 TRUE ~ "None"))
  
  vis_grid = as.data.frame(sim$vis_raster, xy=T) %>%
    select(X = x, Y = y)
  vis_grid$ento_inoculation_rate = apply(vis_grid, 1, sim$calculate_EIR)
  
  ggplot(vis_grid, aes(x=X, y=Y)) +
    geom_raster(aes(fill=ento_inoculation_rate), interpolate=TRUE) +
    scale_fill_viridis_c(option="inferno", limits=c(0, NA)) +
    labs(fill = "EIR") +
    new_scale("fill") +
    geom_point(aes(fill=p_blood_gametocyte, color=Infection), data=humans, pch=21, size=2, stroke=1) +
    labs(fill = "Gametocyte\nprobability") +
    scale_fill_viridis_c(limits=c(0, 1)) +
    scale_color_manual(values=c("New"="tomato",
                                "Historical"="grey",
                                "None"="steelblue")) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(legend.key.height = unit(10, "pt"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())
}

plot_epicurve = function(sim) {
  ymax = sim$history_infections %>%
    filter(source != "Seed") %>%
    count(t_infection) %>%
    pull(n) %>%
    max()
  sim$history_infections %>%
    mutate(source = fct_inorder(source)) %>%
    ggplot(aes(x = t_infection, fill = source)) +
    geom_bar(width = 0.9) +
    coord_cartesian(xlim = c(0, sim$t), ylim = c(0, max(1, ymax))) +
    labs(x = "Infection time",
         fill = NULL,
         y = NULL) +
    scale_fill_brewer(palette = "Set2", na.value = "grey") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
}


#' Plot spread of mosquitoes over distance
plot_mosquito_migration = function(sim) {
  DX = seq(-sim$max_mosquito_flight_range,
           sim$max_mosquito_flight_range,
           length.out = 100)
  DT = seq(1, sim$max_mosquito_lifespan, length.out=7)
  expand.grid(dx=DX, dt=DT) %>%
    mutate(density = sim$mosquito_migration(dx, 0, dt)) %>%
    ggplot(aes(x=dx, y=density, color=dt, group=dt)) +
    geom_line() +
    scale_color_continuous(breaks = DT) +
    labs(x = "Distance (1D)", color = "Days since\nblood meal")
}


plot_mosquito_infectivity = function(sim) {
  tibble(dt = seq(0, sim$max_mosquito_lifespan, length.out=1000),
         survival = sim$mosquito_survival(dt),
         sporozoites = sim$mosquito_sporogony(dt),
         infectivity = survival*sporozoites) %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name, linetype = name)) +
    geom_line() +
    labs(x = "Days since blood meal", y = "Probability", color = "Attribute")
}

plot_human_infectivity = function(sim) {
  tibble(dt = seq(0, 2*sim$duration_human_infectivity, length.out=1000),
         gametocyte_load = sim$human_infectivity(dt),
         immunity = sim$human_immunity(dt)) %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name)) +
    geom_line() +
    labs(title = "Inoculation load and immunity after infection",
         x = "Days since human infection", y = "Probability", color = "Attribute")
}

plot_human_relapse = function(sim) {
  tibble(dt = sim$human_schedule_relapse(10000)) %>%
    ggplot(aes(x = dt)) +
    geom_density(fill = "steelblue") +
    labs(title = "Relapse density",
         x = "Days since infection", y = "Probability")
}