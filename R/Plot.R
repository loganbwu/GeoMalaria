#' Print simulation
#' 
#' @param x Simulation object
#' @param ... Additional arguments
print.Simulation = function(x, ...) {
  print(
    list(
      "t" = x$t,
      "History" = x$linelist
    )
  )
}

#' Plot initial state of a simulation
#' 
#' Valid at any time period, output should not depend on whether simulation is run
#' 
#' @param sim Simulation object
plot_init = function(sim) {
  # Add visible bindings
  ID <- x <- y <- t_infection <- State <- location_proportions <- NULL
  
  # Assemble plotting data
  mosquito_data = as.data.frame(sim$mosquito_raster, xy=T)
  human_data = sim$humans_expand %>%
    mutate(State = factor(case_when(t_infection == 0 ~ "Importation",
                                    TRUE ~ "Susceptible"),
                          levels = c("Susceptible", "Importation"))) %>%
    arrange(State)
  
  ggplot(mapping = aes(x=x, y=y)) +
    geom_raster(data = mosquito_data, aes(fill=layer)) +
    geom_line(data = human_data, aes(color=State, group = ID), alpha=0.6) +
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

#' Plot current state of a simulation
#' 
#' Plots map with epicurve
#' 
#' @param x Simulation object
#' @param t Optional, plot simulation as it was at this time
#' @param ... Additional arguments
plot.Simulation = function(x, t=NULL, ...) {
  p_state = plot_state(x, t, ...)
  p_epicurve = plot_epicurve(x, t, ...)
  p_state / p_epicurve + plot_layout(heights=c(5,1), guides="collect")
}

#' Plot map of simulation
#' 
#' @param sim Simulation object
#' @param t Optional, plot simulation as it was at this time
#' @param ... Unused
plot_state = function(sim, t=NULL, ...) {
  # Add visible bindings
  ID <- x <- y <- t_infection <- ento_inoculation_rate <- location_proportions <- p_blood_gametocyte <- Infection <- NULL
  dot_args = list(...)
  
  if (is.null(t)) {
    t = sim$t
  } else {
    stopifnot("Can only plot a previous state if 'EIR' is added to `log_options`" = t == sim$t | "EIR" %in% names(sim$log))
    stopifnot("Requested t exceeds simulated clock time" = t <= sim$t)
  }
  
  linelist = sim$linelist %>%
    filter(t_infection <= t) %>%
    rev() %>%
    group_by(ID) %>%
    slice(1) %>%
    ungroup() %>%
    full_join(sim$humans_expand[,c("ID", "x", "y", "location_proportions")], by="ID")
  linelist$source[is.na(linelist$source)] = "None"
  
  if ("EIR" %in% names(sim$log)) {
    # Use the closest (most recent) logged EIR to the requested time t
    eir_grid = sim$EIR[sim$EIR$t == max(sim$EIR$t[sim$EIR$t <= t]),]
  } else {
    # Calculate EIR with current state
    eir_grid = sim$vis_grid
    eir_grid$ento_inoculation_rate = apply(sim$vis_grid, 1, sim$calculate_EIR)
  }
  
  # Use the given eir_max if provided in dot arguments, e.g. for an animation
  if ("eir_max" %in% names(dot_args)) {
    eir_max = dot_args$eir_max
  } else {
    eir_max = max(eir_grid$ento_inoculation_rate)
  }
  
  # Create mosquito raster border
  mosquito_border = sim$mosquito_raster
  mosquito_border[] = ifelse(mosquito_border[] > 0, 1, NA)
  mosquito_border = st_as_sf(raster::rasterToPolygons(mosquito_border, dissolve=T))
  
  
  sources = levels(fct_inorder(sim$linelist$source))
  colors = brewer.pal(max(3, length(sources)), "Set2")[seq_along(sources)]
  names(colors) = sources
  colors = c(colors, "None"="grey")
  
  ggplot(eir_grid, aes(x=x, y=y)) +
    geom_raster(aes(fill=ento_inoculation_rate), interpolate=TRUE) +
    geom_sf(data=mosquito_border, aes(x=NULL, y=NULL), color="white", alpha=0, size=0.5) +
    scale_fill_viridis_c(option="inferno", limits=c(0, eir_max)) +
    labs(fill = "EIR") +
    new_scale("fill") +
    geom_line(data=linelist, aes(color=source, group = ID), alpha=0.6) +
    geom_point(data=linelist, aes(color=source, size=location_proportions, stroke=location_proportions)) +
    labs(fill = "Gametocyte\nprobability") +
    scale_size_area(max_size = 2, guide = "none") +
    scale_fill_viridis_c(limits=c(0, 1)) +
    scale_color_manual(values = colors, na.value = "grey") +
    coord_sf() +
    labs(x = NULL, y = NULL, color=NULL) +
    theme_minimal() +
    theme(legend.key.height = unit(10, "pt"),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank())
}

#' Plot epidemic curve
#' 
#' @param sim Simulation object
#' @param t Optional, plot epicurve as it was at this time
#' @param ... Unused
plot_epicurve = function(sim, t=NULL, ...) {
  # Add visible bindings
  t_infection <- NULL
  dot_args = list(...)
  
  if (is.null(t)) {
    t = sim$t
  }
  show_time = t
  
  # Handle dot args
  if ("t_range" %in% names(dot_args)) {
    t_range = dot_args$t_range
  } else {
    t_range = list(min = 0, max = t)
  }
  
  if ("y_max" %in% names(dot_args)) {
    y_max = dot_args$y_max
  } else {
    y_max = with(sim$linelist[sim$linelist$t_infection <= t &
                                  sim$linelist$source != "Seed",],
                   unname(rev(sort(table(t_infection)))[1]))
  }
  
  sources = levels(fct_inorder(sim$linelist$source))
  colors = brewer.pal(max(3, length(sources)), "Set2")[seq_along(sources)]
  names(colors) = sources
  
  sim$linelist %>%
    mutate(source = fct_inorder(source)) %>%
    filter(t_infection <= t) %>%
    ggplot(aes(x = t_infection, fill = source)) +
    geom_vline(xintercept = show_time, color="grey") +
    geom_bar(width = 0.9 * sim$min_dt) +
    coord_cartesian(xlim = c(t_range$min, t_range$max), ylim = c(0, max(1, y_max))) +
    labs(x = "Infection time",
         fill = NULL,
         y = NULL) +
    scale_fill_manual(values = colors, na.value = "grey") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
}

#' Animates a simulation
#' 
#' Plots map with epicurve. Pre-calculates ranges to ensure frames are consistent.
#' 
#' @param sim Simulation object
#' @param t Optional, vector of times to plot
#' @param file Optional, write video to a destination
#' @param ... Additional arguments
plot_anim = function(sim, t=NULL, file=NULL, ...) {
  if (is.null(t)) {
    t_range = list(min = 0,
                   max = sim$t)
    t = seq(t_range$min, t_range$max)
  } else {
    t_range = list(min = min(t),
                   max = max(t))
  }
  
  y_max = with(sim$linelist[sim$linelist$t_infection >= t_range$min &
                                sim$linelist$t_infection <= t_range$max &
                                sim$linelist$source != "Seed",],
                 unname(rev(sort(table(t_infection)))[1]))
  
  eir_max = with(sim$EIR[sim$EIR$t >= t_range$min & sim$EIR$t <= t_range$max,],
                 max(0, max(ento_inoculation_rate)))
  makeplot = function() {
    frames = lapply(t, function(tt) {
      p = plot(sim, tt, t_range=t_range, eir_max=eir_max, y_max=y_max)
      print(p)
    })
  }
  av::av_capture_graphics(makeplot(), file, 1280, 720, res = 144, framerate=15)
}


#' Plot spread of mosquitoes over distance
#' 
#' @param sim Simulation object
plot_mosquito_migration = function(sim) {
  dx = seq(-sim$max_mosquito_flight_range,
           sim$max_mosquito_flight_range,
           length.out = 100)
  dt = seq(1, sim$max_mosquito_lifespan, length.out=7)
  expand.grid(dx=dx, dt=dt) %>%
    mutate(density = sim$mosquito_migration(dx, 0, dt)) %>%
    ggplot(aes(x=dx, y=density, color=dt, group=dt)) +
    geom_line() +
    scale_color_binned(breaks = dt) +
    labs(title = "Distance travelled by mosquitoes",
         x = "Distance (1D shown; actual diffusion is 2D)", y = "Density",
         color = "Days since\nblood meal")
}

#' Plot mosquito infectivity curve
#' 
#' @param sim Simulation object
plot_mosquito_infectivity = function(sim) {
  # Add visible bindings
  survival <- sporozoites <- value <- name <- NULL
  
  tibble(dt = seq(0, sim$max_mosquito_lifespan, length.out=1000),
         survival = sim$mosquito_survival(dt),
         sporozoites = sim$mosquito_sporogony(dt),
         infectivity = survival*sporozoites) %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name, linetype = name)) +
    geom_line() +
    labs(title = "Inoculation probability after blood meal",
         x = "Days since blood meal", y = "Probability", color = "Attribute")
}

#' Plot human infectivity curve
#' 
#' @param sim Simulation object
plot_human_infectivity = function(sim) {
  # Add visible bindings
  value <- name <- NULL
  
  tibble(dt = seq(0, 2*sim$duration_human_infectivity, length.out=1000),
         gametocyte_load = sim$human_infectivity(dt=dt),
         immunity = sim$human_immunity(dt=dt)) %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name)) +
    geom_line() +
    labs(title = "Inoculation load and immunity after infection",
         x = "Days since human infection", y = "Probability", color = "Attribute")
}

#' Plot human relapse risk
#' 
#' @param sim Simulation object
plot_human_relapse = function(sim) {
  tibble(dt = sim$human_schedule_relapse(10000)) %>%
    ggplot(aes(x = dt)) +
    geom_density(fill = "steelblue", alpha = 0.6) +
    scale_x_continuous(limits = c(0, NA)) +
    labs(title = "Relapse density",
         x = "Days since infection", y = "Probability")
}
