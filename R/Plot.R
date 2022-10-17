#' Print simulation
#' 
#' @param x Simulation object
#' @param ... Additional arguments
print.Simulation = function(x, ...) {
  print(
    list(
      "t" = x$t,
      "linelist" = x$linelist,
      "humans" = x$humans,
      "attack_rate" = length(unique(x$linelist$ID)) / nrow(x$humans)
    )
  )
}

#' Plot initial state of a simulation
#' 
#' Valid at any time period, output should not depend on whether simulation is run
#' 
#' @param sim Simulation object
plot_init = function(sim) {
  plot.Simulation(sim, t=0, init=TRUE)
}

#' Plot current state of a simulation
#' 
#' Plots map with epicurve
#' 
#' @param x Simulation object
#' @param t Optional, plot simulation as it was at this time
#' @param init Logical, show initialisation state with mosquito raster
#' @param ... Additional arguments
plot.Simulation = function(x, t=NULL, init=FALSE, ...) {
  if (init) {
    plot_state(x, t, background="mosquito", ...)
  } else {
    p_state = plot_state(x, t, ...)
    p_epicurve = plot_epicurve(x, t, ...) + guides(fill = "none")
    p_state / p_epicurve + plot_layout(heights=c(5,1), guides="collect")
  }
}

#' Plot map of simulation
#' 
#' @param sim Simulation object
#' @param t Optional, plot simulation as it was at this time
#' @param background Either 'EIR' or 'mosquito' 
#' @param ... Unused
plot_state = function(sim, t=NULL, background="EIR", ...) {
  # Add visible bindings
  ID <- x <- y <- t_infection <- location_proportions <- intensity <- NULL
  stopifnot("Invalid argument specified for plot background" = background %in% c("EIR", "mosquito"))
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
  linelist$source = fct_inorder(linelist$source)
  
  # Either use the mosquito raster or EIR calculation as the background raster
  if (background == "mosquito") {
    background = as.data.frame(sim$mosquito_raster, xy=T)
    background_max = max(background$layer)
    background_label = "Mos/km^2"
  } else {
    
    if ("EIR" %in% names(sim$log)) {
      # Use the closest (most recent) logged EIR to the requested time t
      background = sim$EIR[sim$EIR$t == max(sim$EIR$t[sim$EIR$t <= t]), c("x", "y", "ento_inoculation_rate")]
      background_max = max(background$ento_inoculation_rate)
    } else {
      # Calculate EIR with current state
      background = sim$vis_grid
      background$intensity = apply(background, 1, sim$calculate_EIR)
      background_max = max(background$intensity)
    }
    
    # Use the given eir_max if provided in dot arguments, e.g. for an animation
    if ("eir_max" %in% names(dot_args)) {
      background_max = dot_args$eir_max
    }
    background_label = "EIR"
  }
  # Rename raster value because it doesn't matter where it comes from
  names(background)[3] = "intensity"
  
  # Create mosquito raster border
  mosquito_border = sim$mosquito_raster
  mosquito_border[] = ifelse(mosquito_border[] > 0, 1, NA)
  mosquito_border = st_as_sf(raster::rasterToPolygons(mosquito_border, dissolve=T))
  
  sources = levels(fct_inorder(sim$linelist$source))
  colors = brewer.pal(max(3, length(sources)), "Set2")[seq_along(sources)]
  names(colors) = sources
  colors = c(colors, "None"="grey")
  ggplot(background, aes(x=x, y=y)) +
    geom_raster(aes(fill=intensity), interpolate=TRUE) +
    scale_fill_viridis_c(option="inferno", limits=c(0, background_max)) +
    geom_sf(data=mosquito_border, aes(x=NULL, y=NULL), color="grey", fill="transparent") +
    geom_line(data=linelist %>% group_by(ID) %>% filter(n()>1), aes(color=source, group = ID), alpha=0.6) +
    geom_point(data=linelist, aes(color=source, size=location_proportions), alpha=0.6) +
    scale_size_area(max_size = 2, breaks = 0.25*1:4) +
    scale_color_manual(values = colors, na.value = "grey") +
    coord_sf() +
    labs(x=NULL, y=NULL, color="Infection", fill=background_label, size="Proportion\ntime spent") +
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
#' @param t Optional, plot frames for T in [min(t), max(t)
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
                 max(0, max(c(-Inf, ento_inoculation_rate))))
  # makeplot = function() {
  #   frames = lapply(t, function(tt) {
  #     p = plot(sim, tt, t_range=t_range, eir_max=eir_max, y_max=y_max)
  #     print(p)
  #   })
  # }
  # t = head(t, 10)
  temp_dir = tempdir()
  frame_files = file.path(temp_dir, paste0("frame_", t, ".png"))
  if (is.null(file)) {
    file = file.path(temp_dir, "animation.mp4")
    print(paste("Saving animation to", file))
  }
  
  tryCatch({
    # No parallelisation
    pb = progress::progress_bar$new(total = length(t))
    for (i in seq_len(length(t))) {
      pb$tick()
      plot(sim, t[i], t_range=t_range, eir_max=eir_max, y_max=y_max)
      ggsave(frame_files[i], width=1920, height=1080, units="px", dpi=150)
    }
    
    # plot_func = function(i) {
    #   print(paste("Frame", i, "of", length(t)))
    #   plot(sim, t[i], t_range=t_range, eir_max=eir_max, y_max=y_max)
    #   ggsave(frame_files[i], width=1920, height=1080, units="px", dpi=150)
    # }
    # mclapply(seq_len(t), plot_func, mc.cores=8)
    
    # Foreach
    # foreach(
    #   i = seq_along(t), 
    #   .combine = 'c'
    # ) %dopar% {
    #   print(paste("Frame", i, "of", length(t)))
    #   plot(sim, t[i], t_range=t_range, eir_max=eir_max, y_max=y_max)
    #   ggsave(frame_files[i], width=1920, height=1080, units="px")
    # }
    av::av_encode_video(frame_files, file, framerate = 15, verbose = F)
  },
  finally = {
    # print(paste("Removing temporary frames from", temp_dir))
    unlink(frame_files)
  }
  )
  # frames = future_lapply(t, function(tt) {
  #   filename = paste0("~/Documents/GeoMalaria/img/frame_", tt, ".png")
  #   print(filename)
  #   plot(sim, tt, t_range=t_range, eir_max=eir_max, y_max=y_max)
  #   ggsave(filename)
  # })
  # av::av_capture_graphics(future_lapply(t, function(tt) {
  #   p = plot(sim, tt, t_range=t_range, eir_max=eir_max, y_max=y_max)
  #   print(p)
  # }, future.seed=TRUE), file, 1280, 720, res = 144, framerate=15)
  
  # Return the output path
  invisible(file)
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
