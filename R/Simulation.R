library(R6)

#' Malaria simulation class
#' 
#' Agent-based simulation for P. vivax malaria in a spatial context
Simulation = R6Class(
  "Simulation",
  
  #' @field t Current simulation time in days
  #' @field min_dt Smallest time step used so far
  #' @field mosquito_raster Raster object of mosquito counts
  #' @field max_mosquito_lifespan Constrain mosquitoes to an absolute maximum lifespan
  #' @field duration_human_infectivity Constrain infectious period
  #' @field max_mosquito_flight_range Constrain mosquito diffusion distance in km
  #' @field bite_rate Bites per mosquito per density
  #' @field mosquito_death_rate Proportion of mosquitoes that die per day, simulated continuously
  #' @field sporozoite_infection_rate Proportion of infectious mosquito bites that can result in infection
  #' @field p_relapse Probability of an infection scheduling a relapse
  #' @field mean_recovery Recovery parameter mean=1/rate. Could be changed for other distributions or functions
  #' @field vis_raster Raster to visualise any spatial fields on
  #' @field humans Data frame of humans
  #' @field locations Data frame of human locations
  #' @field mosquito_infections Data frame of current mosquito infection events
  #' @field mosquito_travel Constant of mosquito travel in variance of location after one day
  #' @field human_relapse_shape Shape constant of the human relapse distribution
  #' @field human_relapse_rate Rate constant of the human relapse distribution
  #' @field log Character vector of logging options. Options include "linelist", "compartment", and/or "EIR"
  public = list(
    # Constants
    t = 0,
    min_dt = 1,
    mosquito_raster = NULL,
    max_mosquito_lifespan = 0, # calculated on initialisation
    duration_human_infectivity = NULL,
    max_mosquito_flight_range = NULL,
    bite_rate = NULL,
    mosquito_death_rate = NULL,
    sporozoite_infection_rate = 0.75,
    p_relapse = NULL,
    mean_recovery = NULL,
    vis_raster = NULL,
    
    # States
    humans = NULL,
    locations = NULL,
    mosquito_infections = tibble::tibble(x = numeric(),
                                         y = numeric(),
                                         t_inoculation = numeric(),
                                         count = numeric(),
                                         density = numeric()),
    log = list(),
    
    #' @description
    #' Initialize a simulation instance
    #' 
    #' @param humans dataframe of residents in the area
    #' @param mosquito_raster raster of the number of mosquitoes per cell
    #' @param duration_human_infectivity maximum infectious period of a human in days
    #' @param mosquito_death_rate proportion of mosquitoes that die per day
    #' @param bite_rate probability of each mosquito biting a human in a cell
    #' @param p_relapse Proportion of infections that trigger relapses
    #' @param mean_recovery Unused
    #' @param log_options Logging options
    initialize = function(humans,
                          locations,
                          mosquito_raster,
                          duration_human_infectivity,
                          mosquito_death_rate,
                          bite_rate,
                          p_relapse = 0,
                          mean_recovery = 14,
                          log_options = c("linelist", "compartment")) {
      
      # Raster for fine visualisation of continuous fields
      bounds = extent(mosquito_raster)
      self$vis_raster = raster(vals=0, nrows=100, ncols=100,
                               xmn=bounds@xmin, xmx=bounds@xmax,
                               ymn=bounds@ymin, ymx=bounds@ymax)
      
      # Parameters
      self$humans = humans
      # Add columns that allow relapses and recoveries to be scheduled
      if (!"t_relapse" %in% names(self$humans)) {
        self$humans$t_relapse = Inf
      }
      if (!"t_recovery" %in% names(self$humans)) {
        self$humans$t_recovery = Inf
      }
      if (!"p_blood_gametocyte" %in% names(self$humans)) {
        self$humans$p_blood_gametocyte = 0
      }
      
      self$locations = tibble::tibble(locations, gametocyte_load = NA_real_, EIR = NA_real_)
      self$set_mosquito_raster(mosquito_raster)
      self$duration_human_infectivity = duration_human_infectivity
      self$bite_rate = bite_rate
      self$mosquito_death_rate = mosquito_death_rate
      self$p_relapse = p_relapse
      self$mean_recovery = mean_recovery
      
      # Calculate reasonable bounds
      # capture the vast majority of the mosquito lifespan
      while (self$mosquito_survival(self$max_mosquito_lifespan) > 0.01) {
        self$max_mosquito_lifespan = self$max_mosquito_lifespan + 1
      }
      # Capture 95% of the distance travelled at the max lifespan
      self$max_mosquito_flight_range = qnorm(0.975, sd=sqrt(self$mosquito_travel * self$max_mosquito_lifespan))
      
      # Initialise log
      if ("linelist" %in% log_options) {
        self$log$linelist = as_tibble(cbind(
          self$humans[!is.na(self$humans$t_infection), c("ID", "t_infection")],
          source = "Seed"))
      }
      if ("compartment" %in% log_options) {
        self$log$compartment = tibble(
          t = self$t,
          infected = sum(!is.na(self$humans$t_infection)))
      }
      if ("EIR" %in% log_options) {
        self$log$EIR = tibble(
          t = self$t,
          self$vis_grid,
          ento_inoculation_rate = apply(self$vis_grid, 1, self$calculate_EIR)
        )
      }
      
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
      self$min_dt = min(self$min_dt, dt)
      
      ## Calculate current state of human infectivity and immunity
      # If a human is due to recover, wipe the infection
      recovery_IDs = which(self$humans$t_recovery <= self$t)
      self$humans$t_infection[recovery_IDs] = NA # wipe infection
      self$humans$t_recovery[recovery_IDs] = Inf # unschedule recovery
      # Otherwise calculate their disease state
      self$humans$p_blood_gametocyte = self$human_infectivity(self$humans, self$t)
      self$humans$immunity = self$human_immunity(self$humans, self$t)
      
      # Sum total human gametocyte load per location
      self$locations$gametocyte_load = 0
      infected = self$humans[self$humans$p_blood_gametocyte > 0,]
      for (i in seq_len(nrow(infected))) {
        with(infected, {
          ix = location_ix[[i]]
          self$locations$gametocyte_load[ix] = self$locations$gametocyte_load[ix] +
            p_blood_gametocyte[i] * # adjust by human infectivity
            location_proportions[[i]] # weight by human time spent
        })
      }
      
      ## Calculate state of infected mosquitoes
      # Remove expired events
      self$mosquito_infections = self$mosquito_infections[
        self$t - self$mosquito_infections$t_inoculation <= self$max_mosquito_lifespan,]
      # Calculate current infectivity of infected mosquito clouds integrated over space
      # i.e., number of live mosquitoes with mature sporozoites from this event
      self$mosquito_infections$density = self$mosquito_infections$count *
        self$mosquito_survival(self$t - self$mosquito_infections$t_inoculation) *
        self$mosquito_sporogony(self$t - self$mosquito_infections$t_inoculation)
      
      # Expose locations
      self$locations$EIR = apply(self$locations, 1, self$calculate_EIR)
      
      # Expose and infect humans who are not immune or are due to relapse
      susceptible = self$humans[self$humans$immunity < 1 | self$humans$t_relapse <= self$t,]
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
        susceptible_EIR = numeric() # otherwise it'll be undefined
      }
      # shortcut for P(x > 0), x ~Pois(EIR*dt) OR relapse
      infect_ix = runif(nrow(susceptible)) > exp(-susceptible_EIR * dt) |
        susceptible$t_relapse <= self$t
      is_relapse = (susceptible$t_relapse <= self$t)[infect_ix]
      infect_IDs = susceptible$ID[infect_ix]
      
      ## Update human population
      self$humans$t_infection[infect_IDs] = self$t
      # clear and reschedule relapses and recoveries
      self$humans$t_relapse[infect_IDs] = self$t + self$human_schedule_relapse(length(infect_IDs))
      self$humans$t_recovery[infect_IDs] = self$t + self$human_schedule_recovery(length(infect_IDs))
      
      ## Update mosquito population
      # Expose mosquitoes at locations that have both gametocytes and mosquitoes
      new_mosquito_infections = self$locations[self$locations$gametocyte_load > 0 &
                                                 self$locations$mosquito_density > 0,
                                               c("x", "y", "gametocyte_load", "mosquito_density")]
      if (nrow(new_mosquito_infections) > 0) {
        # Infect mosquitoes at rate (assume no depletion of the susceptible compartment)
        new_mosquito_infections$count = new_mosquito_infections$gametocyte_load *
          new_mosquito_infections$mosquito_density * # N.B. local density could potentially be recomputed
          self$bite_rate *
          dt
        new_mosquito_infections$density = 0 # gets overwritten later anyway but prevent NAs
        new_mosquito_infections$t_inoculation = self$t
        new_mosquito_infections$gametocyte_load <- new_mosquito_infections$mosquito_density <- NULL
        # Add new potential mosquitoes to list
        self$mosquito_infections = bind_rows(self$mosquito_infections, new_mosquito_infections)
      }
      
      # Log observations
      if ("linelist" %in% names(self$log)) {
        self$log$linelist = bind_rows(
          self$log$linelist,
          tibble::tibble(ID = infect_IDs,
                         t_infection = self$t,
                         source = as.character(ifelse(is_relapse, "Relapse", "Transmission"))))
      }
      if ("compartment" %in% names(self$log)) {
        self$log$compartment = bind_rows(
          self$log$compartment,
          tibble(t = self$t,
                 infected = sum(!is.na(self$humans$t_infection))))
      }
      if ("EIR" %in% names(self$log)) {
        self$log$EIR = bind_rows(
          self$log$EIR,
          tibble(
            t = self$t,
            self$vis_grid,
            ento_inoculation_rate = apply(self$vis_grid, 1, self$calculate_EIR)
          )
        )
      }
      
      invisible(self)
    },
    
    #' Recompute location mosquito densities from a raster
    #' 
    #' Updates simulation's raster and densities at each location
    #' 
    #' @param mosquito_raster Raster of mosquito densities
    set_mosquito_raster = function(mosquito_raster) {
      self$locations$mosquito_density = my_extract(mosquito_raster, self$locations[c("x", "y")])
      self$mosquito_raster = mosquito_raster
      invisible(mosquito_raster)
    },
    
    #' Entomological inoculation rate
    #' 
    #' Calculate EIR at a location with x and y elements
    #' 
    #' @param loc List or dataframe with `x` and `y` columns
    #' 
    #' @return Numeric vector of EIR
    calculate_EIR = function(loc) {
      # Calculate total EIR=HBR*SIR from all mosquito infection events
      with(
        self$mosquito_infections,
        sum(density *
              self$bite_rate *
              # Disperse load over total cloud area (TODO: check integral=1)
              self$mosquito_migration(x-loc["x"], y-loc["y"], self$t-t_inoculation) *
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
    mosquito_migration = function(dx, dy, dt) {
      exp(dnorm(dx, sd=sqrt(self$mosquito_travel * dt), log=T) + 
            dnorm(dy, sd=sqrt(self$mosquito_travel * dt), log=T))
    },
    
    #' Lifecycle stages for P. vivax in humans after inoculation
    #' 
    #' Calculates a vector of probability of infecting mosquito with gametocytes in the grid cell.
    #' Specify either `humans` AND `t`, OR `dt`
    #' 
    #' @param humans Dataframe of humans with a `t_infection` column
    #' @param t Simulation time
    #' @param dt Time since infection
    human_infectivity = function(humans=NULL, t=NULL, dt=NULL) {
      stopifnot("Function needs (humans AND t) OR dt" = xor(!is.null(humans) & !is.null(t), !is.null(dt)))
      if (is.null(dt)) {
        dt = t - humans$t_infection
      }
      
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
    #' after decay, everyone retains partial protection.
    #' 
    #' Specify either `humans` AND `t`, OR `dt`
    #' 
    #' @param humans Dataframe of humans with a `t_infection` column
    #' @param t Simulation time
    #' @param dt Time since infection
    human_immunity = function(humans=NULL, t=NULL, dt=NULL) {
      stopifnot("Function needs (humans AND t) OR dt" = xor(!is.null(humans) & !is.null(t), !is.null(dt)))
      if (is.null(dt)) {
        dt = t - humans$t_infection
      }
      
      immunity = approx(c(0, self$duration_human_infectivity, 2*self$duration_human_infectivity),
                        c(1, 1, 0),
                        dt,
                        yleft = 0, yright = 1)$y
      replace_na(immunity, 0)
    },
    
    
    human_relapse_shape = 1,
    human_relapse_rate = 1/30,
    #' Human relapse distribution
    #' 
    #' Schedule a relapse infection for number of days after the primary infection.
    #' Returns `Inf` value if no relapse.
    #' shape = (mean / sd)^2
    #' rate =  mean / sd^2
    #' 
    #' @param n Number of people to schedule relapses for
    human_schedule_relapse = function(n) {
      ifelse(runif(n) < self$p_relapse,
             14 + rgamma(n, self$human_relapse_shape, self$human_relapse_rate),
             Inf)
    },
    
    #' Human recovery distribution
    #' 
    #' Schedule a recovery for number of days after the primary infection. Can
    #' also be used to schedule the success of treatment.
    #' Returns `Inf` value if no recovery
    #' 
    #' @param n Number of people to schedule relapses for
    human_schedule_recovery = function(n) {
      # rexp(n, rate = 1 / self$mean_recovery)
      Inf
    }
  ),
  
  #' @field humans_collapse Represent humans with only one set of coordinates
  #' @field humans_expand Represent humans at all of their coordinates
  #' @field linelist Safely access linelist if logged
  #' @field compartment Safely access compartments if logged
  #' @field EIR Entomological inoculation rate by coordinates over time
  #' @field vis_grid Data frame version of `vis_raster`
  active = list(
    
    linelist = function() {
      stopifnot("Tried to access linelist logs but 'linelist' not specified in `log_options`" = "linelist" %in% names(self$log))
      self$log$linelist
    },
    
    compartment = function() {
      stopifnot("Tried to access compartment size logs but 'compartment' not specified in `log_options`" = "compartment" %in% names(self$log))
      self$log$compartment
    },
    
    EIR = function() {
      stopifnot("Tried to access EIR logs but 'EIR' not specified in `log_options`" = "EIR" %in% names(self$log))
      self$log$EIR
    },
    
    vis_grid = function() {
      as.data.frame(self$vis_raster, xy=T)[, c("x", "y")]
    },
    
    #' Representation of humans with only one set of coordinates
    #' 
    #' Just places people at their highest proportion location
    humans_collapse = function() {
      self$humans %>%
        mutate(location_ix = unlist(map(location_ix, first)),
               x = self$locations$x[location_ix],
               y = self$locations$y[location_ix])
    },
    
    #' Represent humans at all of their coordinates
    humans_expand = function() {
      self$humans %>%
        unnest(cols=c(location_ix, location_proportions)) %>%
        mutate(x = self$locations$x[location_ix],
               y = self$locations$y[location_ix])
    }
  )
)
