library(R6)

#' Malaria simulation class
#' 
#' Agent-based simulation for P. vivax malaria in a spatial context
Simulation = R6Class(
  "Simulation",
  
  #' @field t Current simulation time in days
  #' @field mosquito_raster Raster object of mosquito counts
  #' @field max_mosquito_lifespan Constrain mosquitoes to an absolute maximum lifespan
  #' @field duration_human_infectivity Constrain infectious period
  #' @field max_mosquito_flight_range Constrain mosquito diffusion distance in km
  #' @field bite_rate Bites per mosquito per density
  public = list(
    # Constants
    t = 0,
    mosquito_raster = NULL,
    max_mosquito_lifespan = 0,
    duration_human_infectivity = NULL,
    max_mosquito_flight_range = NULL,
    bite_rate = NULL,
    mosquito_death_rate = NULL,
    sporozoite_infection_rate = 0.75,
    p_relapse = NULL,
    
    # States
    humans = NULL,
    locations = NULL,
    vis_raster = NULL,
    mosquito_infections = tibble::tibble(X = numeric(),
                                 Y = numeric(),
                                 t_inoculation = numeric()),
    humans_infections = tibble::tibble(),
    
    # Outbreak history
    history_infections = NULL,
    
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
                          bite_rate,
                          p_relapse = 0) { # proportion of infections that relapse including relapses themselves
      
      # Raster for fine visualisation of continuous fields
      bounds = extent(mosquito_raster)
      self$vis_raster = raster(vals=0, nrows=100, ncols=100,
                                  xmn=bounds@xmin, xmx=bounds@xmax,
                                  ymn=bounds@ymin, ymx=bounds@ymax)
      
      # Parameters
      self$humans = humans
      if (!"t_relapse" %in% names(self$humans)) {
        self$humans$t_relapse = Inf
      }
      self$locations = tibble::tibble(locations, gametocyte_load = NA_real_, EIR = NA_real_)
      self$mosquito_infections = tibble::tibble(X = numeric(), Y = numeric(), infected_count = integer(), t_inoculation = numeric())
      self$set_mosquito_raster(mosquito_raster)
      self$duration_human_infectivity = duration_human_infectivity
      self$bite_rate = bite_rate
      self$mosquito_death_rate = mosquito_death_rate
      self$p_relapse = p_relapse
      
      # Calculate reasonable mosquito bounds
      # capture the vast majority of the mosquito lifespan
      while (self$mosquito_survival(self$max_mosquito_lifespan) > 0.01) {
        self$max_mosquito_lifespan = self$max_mosquito_lifespan + 1
      }
      # Capture 95% of the distance travelled at the max lifespan
      self$max_mosquito_flight_range = qnorm(0.975, sd=sqrt(self$mosquito_travel * self$max_mosquito_lifespan))
      
      # Initialise history
      self$history_infections = self$humans %>%
        filter(!is.na(t_infection)) %>%
        select(ID, t_infection) %>%
        mutate(source = "Seed")
      
      invisible(self)
    },
    
    
    #' Iterate simulation
    #' 
    #' Advance simulation by `dt`. Calculates the current state of humans, infects
    #' mosquitoes and humans simultaneously, updates humans and mosquitoes, then
    #' stores observations for analysis.
    #' 
    #' @param dt Time to advance simulation by in days
    #' @param debug Optional, print debug information
    iterate = function(dt, debug=FALSE) {
      self$t = self$t + dt
      
      # Calculate current state of human infectivity and immunity
      self$humans = self$humans %>%
        mutate(p_blood_gametocyte = self$human_infectivity(self$t - t_infection),
               immunity = self$human_immunity(self$t - t_infection))
      
      # Calculate total human gametocyte load per location
      self$locations$gametocyte_load = 0
      infected = self$humans %>%
        filter(p_blood_gametocyte > 0)
      for (i in seq_len(nrow(infected))) {
        with(infected, {
          ix = location_ix[[i]]
          self$locations$gametocyte_load[ix] = self$locations$gametocyte_load[ix] +
            p_blood_gametocyte[i] * # adjust by human infectivity
            location_proportions[[i]] # weight by human time spent
        })
      }
      
      # Calculate state of infected mosquitoes
      self$mosquito_infections = self$mosquito_infections %>%
        # Remove expired events
        filter(self$t - t_inoculation <= self$max_mosquito_lifespan) %>%
        # Calculate current infectivity of infected mosquito clouds integrated over space   # Number of live mosquitoes with mature sporozoites from this event
        mutate(infectious_count = infected_count *
                 self$mosquito_survival(self$t - t_inoculation) *
                 self$mosquito_sporogony(self$t - t_inoculation))
      
      # Expose locations
      self$locations$EIR = apply(self$locations, 1, self$calculate_EIR)
      
      # Expose and infect humans
      susceptible = self$humans %>%
        # People who are not immune or are due to relapse
        filter(immunity < 1 | t_relapse <= self$t)
      if (nrow(susceptible) > 0) {
        susceptible_EIR = apply(
          susceptible,
          1,
          function(human) {
            ix = human$location_ix
            # get weighted average of EIR at the person's locations
            sum(self$locations$EIR[ix] * human$location_proportions)
          }
        ) * (1 - susceptible$immunity)
      } else {
        susceptible_EIR = numeric()
      }
      # shortcut for P(x > 0), x ~Pois(EIR*dt) OR relapse
      infect_ix = runif(nrow(susceptible)) > exp(-susceptible_EIR * dt) |
        susceptible$t_relapse <= self$t
      is_relapse = (susceptible$t_relapse <= self$t)[infect_ix]
      infect_IDs = susceptible$ID[infect_ix]
      
      # Update human population
      self$humans$t_infection[infect_IDs] = self$t
      # clear and reschedule relapses
      self$humans$t_relapse[infect_IDs] = self$t + self$human_schedule_relapse(length(infect_IDs))
      
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
      
      # Take observations
      self$history_infections = bind_rows(
        self$history_infections,
        tibble::tibble(ID = infect_IDs,
               t_infection = self$t,
               source = as.character(ifelse(is_relapse, "Relapse", "Transmission")))
      )
      
      invisible(self)
    },
    
    #' Recompute location mosquito densities from a raster
    #' 
    #' Updates simulation's raster and densities at each location
    #' 
    #' @param mosquito_raster Raster of mosquito densities
    set_mosquito_raster = function(mosquito_raster) {
      self$locations$mosquito_density = my_extract(mosquito_raster, self$locations[c("X", "Y")])
      self$mosquito_raster = mosquito_raster
      invisible(mosquito_raster)
    },
    
    #' Entomological inoculation rate
    #' 
    #' Calculate EIR at a location with X and Y elements
    #' 
    #' @param loc List or dataframe with `X` and `Y` columns
    #' 
    #' @return Numeric vector of EIR
    calculate_EIR = function(loc) {
      # Calculate total EIR=HBR*SIR from all mosquito infection events
      with(
        self$mosquito_infections,
        sum(infectious_count *
              self$bite_rate *
              # Disperse load over total cloud area (TODO: check integral=1)
              self$mosquito_migration(X-loc["X"], Y-loc["Y"], self$t-t_inoculation) *
              self$sporozoite_infection_rate)
      )
    },
    
    #' Proportion of mosquitoes that survive over time
    #' 
    #' @param dt Time since mosquito inoculation in days
    mosquito_survival = function(dt) {
      survival = dexp(dt, rate=self$mosquito_death_rate) / self$mosquito_death_rate
    },
    
    #' Proportion of mosquitoes with sporozoites over time
    #' 
    #' @param dt Time since mosquito inoculation in days
    mosquito_sporogony = function(dt) {
      approx(c(8, 10), c(0, 1), dt, yleft=0, yright=1)$y
    },
    
    mosquito_travel = 1^2,
    #' Mosquito density
    #' 
    #' Mosquito density calculated if they started at a point and traveled `dx` and
    #' `dy` in time `dt`. Variance is constant at (km/day)^2. Uses the Laplace isometric
    #' diffusion solution.
    #' 
    #' @param dx Distance traveled horizontally in km
    #' @param dy Distance traveled vertically
    #' @param dt Time since origin in days
    #' 
    #' @return Relative mosquito density
    mosquito_migration = function(dx, dy, dt) {
      exp(dnorm(dx, sd=sqrt(self$mosquito_travel * dt), log=T) + 
            dnorm(dy, sd=sqrt(self$mosquito_travel * dt), log=T))
    },
    
    #' Lifecycle stages for P. vivax in humans after inoculation
    #'
    #' @param dt Time since inoculation in days
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
    #' 
    #' @param dt Time since infection in days
    human_immunity = function(dt) {
      immunity = approx(c(0, self$duration_human_infectivity, 2*self$duration_human_infectivity),
                        c(1, 1, 1),
                        dt,
                        yleft = 0, yright = 1)$y
      replace_na(immunity, 0)
    },
    
    
    human_relapse_shape = 1,
    human_relapse_rate = 1/30,
    #' Human relapse distribution
    #' 
    #' Schedule a relapse infection for number of days after the primary infection.
    #' NA value is no relapse.
    #' shape = (mean / sd)^2
    #' rate =  mean / sd^2
    #' 
    #' @param n Number of people to schedule relapses for
    #' 
    #' @return Time until scheduled relapse, `Inf` if no relapse
    human_schedule_relapse = function(n) {
      ifelse(runif(n) < self$p_relapse,
             14 + rgamma(n, self$human_relapse_shape, self$human_relapse_rate),
             Inf)
    }
  ),
  
  
  active = list(
    
    #' Representation of humans with only one set of coordinates
    #' 
    #' Just places people at their highest proportion location
    #' 
    #' @return Tibble of humans
    humans_collapse = function() {
      self$humans %>%
        mutate(location_ix = unlist(map(location_ix, first)),
               X = self$locations$X[location_ix],
               Y = self$locations$Y[location_ix])
    },
    
    #' Represent humans at all of their coordinates
    #' 
    #' @return Tibble of humans duplicated if at multiple locations
    humans_expand = function() {
      self$humans %>%
        unnest(cols=c(location_ix, location_proportions)) %>%
        mutate(X = self$locations$X[location_ix],
               Y = self$locations$Y[location_ix])
    }
  )
)
