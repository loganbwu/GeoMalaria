library(R6)

Simulation = R6Class(
  "Simulation",
  public = list(
    # Constants
    t = 0,
    mosquito_raster = NULL,
    n_humans = NULL,
    max_mosquito_lifespan = 0,
    duration_human_infectivity = NULL,
    max_mosquito_flight_range = NULL,
    bite_rate = NULL,
    mosquito_death_rate = NULL,
    sporozoite_infection_rate = 0.75,
    
    # States
    humans = NULL,
    locations = NULL,
    mosquito_infections = tibble(X = numeric(),
                                 Y = numeric(),
                                 t_inoculation = numeric()),
    humans_infections = tibble(),
    
    #' Initialize a simulation instance
    #' 
    #' @param humans dataframe of residents in the area
    #' @param mosquito_raster raster of the number of mosquitoes per cell
    #' @param duration_human_infectivity maximum infectious period of a human in days
    #' @param mosquito_death_rate proportion of mosquitoes that die per day
    #' @param bite_rate probability of each mosquito biting a human in a cell
    initialize = function(humans,
                          locations,
                          mosquito_raster,
                          duration_human_infectivity, # days
                          mosquito_death_rate, # per day
                          bite_rate) {
      
      # Raster for fine visualisation of continuous fields
      bounds = extent(mosquito_raster)
      private$vis_raster = raster(vals=0, nrows=100, ncols=100,
                                  xmn=bounds@xmin, xmx=bounds@xmax,
                                  ymn=bounds@ymin, ymx=bounds@ymax, 
                                  crs="+units=km")
      
      # Parameters
      self$humans = humans
      self$locations = tibble(locations, gametocyte_load = NA_real_, EIR = NA_real_)
      self$mosquito_infections = tibble(X = numeric(), Y = numeric(), infected_count = integer(), t_inoculation = numeric())
      self$set_mosquito_raster(mosquito_raster)
      self$duration_human_infectivity = duration_human_infectivity
      self$bite_rate = bite_rate
      self$mosquito_death_rate = mosquito_death_rate
      
      # Calculate reasonable mosquito bounds
      # capture the vast majority of the mosquito lifespan
      while (private$mosquito_survival(self$max_mosquito_lifespan) > 0.01) {
        self$max_mosquito_lifespan = self$max_mosquito_lifespan + 1
      }
      # Capture 95% of the distance travelled at the max lifespan
      self$max_mosquito_flight_range = qnorm(0.975, sd=sqrt(private$mosquito_travel * self$max_mosquito_lifespan))
      
      invisible(self)
    },
    
    
    iterate = function(dt, debug=FALSE) {
      self$t = self$t + dt
      
      # Calculate current state of human infectivity and immunity
      self$humans = self$humans %>%
        mutate(blood_gametocyte_prob = private$human_infectivity(self$t - t_infection),
               immunity = private$human_immunity(self$t - t_infection))
      
      # Calculate total human gametocyte load per location
      self$locations$gametocyte_load = 0
      infected = self$humans %>%
        filter(blood_gametocyte_prob > 0)
      for (i in seq_len(nrow(infected))) {
        with(infected, {
          ix = location_ix[[i]]
          self$locations$gametocyte_load[ix] = self$locations$gametocyte_load[ix] +
            blood_gametocyte_prob[i] * # adjust by human infectivity
            location_proportions[[i]] # weight by human time spent
        })
      }
      
      # Calculate state of infected mosquitoes
      self$mosquito_infections = self$mosquito_infections %>%
        # Remove expired events
        filter(self$t - t_inoculation <= self$max_mosquito_lifespan) %>%
        # Calculate current infectivity of infected mosquito clouds integrated over space   # Number of live mosquitoes with mature sporozoites from this event
        mutate(infectious_count = infected_count *
                 private$mosquito_survival(self$t - t_inoculation) *
                 private$mosquito_sporogony(self$t - t_inoculation))
      
      # Expose locations
      self$locations$EIR = apply(locations, 1, private$calculate_EIR)
      
      # Expose and infect humans
      susceptible = mysim$humans %>%
        filter(immunity < 1)
      susceptible_EIR = apply(
        susceptible,
        1,
        function(human) {
          ix = human$location_ix
          # get weighted average of EIR at the person's locations
          sum(self$locations$EIR[ix] * human$location_proportions)
        }
      ) * (1 - susceptible$immunity)
      # shortcut for P(x > 0), x ~Pois(EIR*dt)
      infect_IDs = susceptible$ID[runif(nrow(susceptible)) > exp(-susceptible_EIR * dt)]
      
      # Update human population
      self$humans$t_infection[infect_IDs] = self$t
      
      # Update mosquito population
      # Expose mosquitoes at human locations
      new_mosquito_infections = self$locations %>%
        filter(gametocyte_load > 0) %>%
        # Infect mosquitoes at rate (assume no susc. depletion)
        # Number of new infected mosquitoes generated by this human at this step
        mutate(infected_count = gametocyte_load *
                 mosquito_density * # N.B. local density could potentially be recomputed
                 self$bite_rate *
                 dt,
               infectious_count = 0, # gets overwritten later anyway but prevent NAs
               t_inoculation = self$t) %>%
        filter(infected_count > 0) %>% # remove zero mosquito/gametocyte events
        select(X, Y, t_inoculation, infected_count, infectious_count)
      # Add potential mosquitoes to list
      self$mosquito_infections = bind_rows(self$mosquito_infections, new_mosquito_infections)
      
      invisible(self)
    },
    
    plot = function() {
      p_state = self$plot_state()
      p_epicurve = self$plot_epicurve()
      p_state / p_epicurve + plot_layout(heights=c(5,1), guides="collect")
    },
    
    
    plot_init = function() {
      mosquito_data = as.data.frame(self$mosquito_raster, xy=T) %>%
        rename(X = x, Y = y)
      humans = self$humans_expand %>%
        mutate(State = case_when(t_infection == 0 ~ "Importation",
                                 TRUE ~ "Susceptible")) %>%
        arrange(t_infection)
      ggplot(mapping = aes(x=X, y=Y)) +
        geom_raster(data = mosquito_data, aes(fill=count)) +
        geom_point(data = humans, aes(color=State, size=location_proportions), alpha=0.6) +
        coord_equal() +
        scale_size_area(max_size = 3) +
        scale_fill_viridis_c(option = "magma") +
        scale_color_manual(values = c("Importation"="tomato",
                                      "Susceptible"="grey")) +
        labs(title = "Initial state",
             fill = "Mos/km^2", size = "Proportion\ntime spent",
             x = NULL, y = NULL) +
        theme_minimal() +
        theme(legend.key.height = unit(10, "pt"),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.background = element_blank(),
              panel.grid = element_blank())
    },
    
    
    plot_state = function() {
      humans = self$humans_collapse %>%
        mutate(Infection = case_when(t_infection == self$t ~ "New",
                                     t_infection < self$t ~ "Historical",
                                     TRUE ~ "None"))
      
      vis_grid = as.data.frame(private$vis_raster, xy=T) %>%
        select(X = x, Y = y)
      print(self$mosquito_infections %>%
              filter(is.na(infectious_count)))
      vis_grid$ento_inoculation_rate = apply(vis_grid, 1, private$calculate_EIR)
      
      ggplot(vis_grid, aes(x=X, y=Y)) +
        geom_raster(aes(fill=ento_inoculation_rate), interpolate=TRUE) +
        scale_fill_viridis_c(option="inferno", limits=c(0, NA)) +
        labs(fill = "EIR") +
        new_scale("fill") +
        geom_point(aes(fill=blood_gametocyte_prob, color=Infection), data=humans, pch=21, size=2, stroke=1) +
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
    },
    
    plot_epicurve = function() {
      self$humans %>%
        filter(t_infection > 0) %>%
        ggplot(aes(x = t_infection)) +
        geom_bar(width = 0.9) +
        coord_cartesian(xlim = c(0, self$t)) +
        labs(x = "Infection time",
             y = NULL) +
        theme_minimal() +
        theme(panel.grid.minor = element_blank())
    },
    
    #' Plot spread of mosquitoes over distance
    plot_mosquito_migration = function() {
      DX = seq(-self$max_mosquito_flight_range,
               self$max_mosquito_flight_range,
               length.out = 100)
      DT = seq(1, self$max_mosquito_lifespan, length.out=7)
      expand.grid(dx=DX, dt=DT) %>%
        mutate(density = private$mosquito_migration(dx, 0, dt)) %>%
        ggplot(aes(x=dx, y=density, color=dt, group=dt)) +
        geom_line() +
        scale_color_continuous(breaks = DT) +
        labs(x = "Distance (1D)", color = "Days since\nblood meal")
    },
    
    
    plot_mosquito_infectivity = function() {
      tibble(dt = seq(0, self$max_mosquito_lifespan, length.out=1000),
             survival = private$mosquito_survival(dt),
             sporozoites = private$mosquito_sporogony(dt),
             infectivity = survival*sporozoites) %>%
        pivot_longer(cols = -dt) %>%
        ggplot(aes(x = dt, y = value, color = name, linetype = name)) +
        geom_line() +
        labs(x = "Days since blood meal", y = "Probability", color = "Attribute")
    },
    
    plot_human_infectivity = function() {
      tibble(dt = seq(0, 2*self$duration_human_infectivity, length.out=1000),
             gametocyte_load = private$human_infectivity(dt),
             immunity = private$human_immunity(dt)) %>%
        pivot_longer(cols = -dt) %>%
        ggplot(aes(x = dt, y = value, color = name)) +
        geom_line() +
        labs(title = "Inoculation load and immunity after infection",
             x = "Days since human infection", y = "Probability", color = "Attribute")
    },
    
    #' Recompute location mosquito densities from a raster
    #' 
    #' @return mosquito_raster for further assignment
    set_mosquito_raster = function(mosquito_raster) {
      self$locations$mosquito_density = my_extract(mosquito_raster, self$locations[c("X", "Y")])
      self$mosquito_raster = mosquito_raster
      invisible(mosquito_raster)
    }
  ),
  
  
  
  private = list(
    env_raster = NULL,
    vis_raster = NULL,
    
    #' Calculate EIR at a location with X and Y elements
    calculate_EIR = function(loc) {
      # Calculate total EIR=HBR*SIR from all mosquito infection events
      with(
        self$mosquito_infections,
        sum(infectious_count *
              self$bite_rate *
              # Disperse load over total cloud area (TODO: check integral=1)
              private$mosquito_migration(X-loc["X"], Y-loc["Y"], self$t-t_inoculation) *
              self$sporozoite_infection_rate)
      )
    },
    
    #' Proportion of mosquitoes that survive over time
    mosquito_survival = function(dt) {
      survival = dexp(dt, rate=self$mosquito_death_rate) / self$mosquito_death_rate
    },
    
    #' Proportion of mosquitoes with sporozoites over time
    mosquito_sporogony = function(dt) {
      approx(c(8, 10), c(0, 1), dt, yleft=0, yright=1)$y
    },
    
    # 1-D variance of mosquito travel in (km/day)^2
    mosquito_travel = 1^2,
    mosquito_migration = function(dx, dy, dt) {
      exp(dnorm(dx, sd=sqrt(private$mosquito_travel * dt), log=T) + 
            dnorm(dy, sd=sqrt(private$mosquito_travel * dt), log=T))
    },
    
    #' Lifecycle stages for P. vivax in humans after inoculation
    #'
    #' @param t days since inoculation
    #' 
    #' @return vector of probability of infecting mosquito with gametocytes in the grid cell
    human_infectivity = function(dt) {
      infectivity = approx(c(12, 17, 20, self$duration_human_infectivity),
                           c(0, 1, 1, 0),
                           dt,
                           yleft = 0, yright = 0)$y
      replace_na(infectivity, 0)
    },
    
    #' Human immunity profile over time
    #' 
    #' Full immunity then decays to 50% 
    #' scales # bites not # susceptible people. E.g. no one retains complete immunity
    #' after decay, everyone retains partial protection
    human_immunity = function(dt) {
      immunity = approx(c(0, self$duration_human_infectivity, 2*self$duration_human_infectivity),
                        c(1, 1, 0.5),
                        dt,
                        yleft = 0, yright = 0.5)$y
      replace_na(immunity, 0)
    }
  ),
  
  
  active = list(
    
    #' Representation of humans with only one set of coordinates
    #' 
    #' Just places people at their highest proportion location
    humans_collapse = function() {
      self$humans %>%
        mutate(location_ix = unlist(map(location_ix, first)),
               X = self$locations$X[location_ix],
               Y = self$locations$Y[location_ix])
    },
    
    #' Represent humans at all of their coordinates
    humans_expand = function() {
      self$humans %>%
        unnest(cols=c(location_ix, location_proportions)) %>%
        mutate(X = self$locations$X[location_ix],
               Y = self$locations$Y[location_ix])
    }
  )
)
