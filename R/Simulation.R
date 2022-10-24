library(R6)

#' Malaria simulation class
#' 
#' Agent-based simulation for P. vivax malaria in a spatial context
Simulation = R6Class(
  "Simulation",
  
  #' @field t Current simulation time in days
  #' @field min_dt Smallest time step used so far
  #' @field mosquito_raster Raster object of mosquito counts
  #' @field sporozoite_infection_rate Proportion of infectious mosquito bites that can result in infection
  #' @field p_relapse Probability of an infection scheduling a relapse
  #' @field mean_recovery Recovery parameter mean=1/rate. Could be changed for other distributions or functions
  #' @field vis_raster Raster to visualise any spatial fields on
  #' @field population Data frame of population
  #' @field locations Data frame of human locations
  #' @field mosquito_infections Data frame of current mosquito infection events
  #' @field mosquito_travel Constant of mosquito travel in variance of location after one day
  #' @field relapse_shape Shape constant of the human relapse distribution
  #' @field relapse_rate Rate constant of the human relapse distribution
  #' @field log Character vector of logging options. Options include "linelist", "compartment", and/or "EIR"
  public = list(
    # Constants
    t = 0,
    min_dt = Inf,
    mosquito_raster = NULL,
    mosquito = NULL,
    human = NULL,
    sporozoite_infection_rate = 0.75,
    p_relapse = NULL,
    mean_recovery = NULL,
    vis_raster = NULL,
    
    # States
    population = NULL,
    locations = NULL,
    mosquito_infections = tibble::tibble(x = numeric(),
                                         y = numeric(),
                                         t_inoculation = numeric(),
                                         count = numeric(),
                                         density = numeric()),
    log = list(),
    
    #' @description Initialize a simulation instance
    #' 
    #' @param population dataframe of residents in the area
    #' @param locations dataframe of human locations
    #' @param mosquito class for mosquitoes
    #' @param mosquito_raster raster of the number of mosquitoes per cell
    #' @param mosquito_death_rate proportion of mosquitoes that die per day
    #' @param bite_rate probability of each mosquito biting a human in a cell
    #' @param p_relapse Proportion of infections that trigger relapses
    #' @param mean_recovery Unused
    #' @param log_options Logging options
    initialize = function(population,
                          locations,
                          mosquito = Mosquito$new(),
                          mosquito_raster,
                          human = Human$new(),
                          p_relapse = 0,
                          mean_recovery = 14,
                          log_options = c("linelist", "compartment")) {
      
      # Raster for fine visualisation of continuous fields
      bounds = extent(mosquito_raster)
      self$vis_raster = raster(vals=0, nrows=100, ncols=100,
                               xmn=bounds@xmin, xmx=bounds@xmax,
                               ymn=bounds@ymin, ymx=bounds@ymax)
      
      # Verify that all locations have defined mosquito amounts
      stopifnot("Locations not covered by raster" = any(
        !is.na(my_extract(mosquito_raster, locations[c("x", "y")]))
      ))
      
      # Parameters
      self$human = human
      self$population = population
      # Add columns for the internal state of population
      if (!"t_relapse" %in% names(self$population)) {
        self$population$t_relapse = Inf
      }
      if (!"t_recovery" %in% names(self$population)) {
        self$population$t_recovery = Inf
      }
      if (!"p_blood_gametocyte" %in% names(self$population)) {
        self$population$p_blood_gametocyte = 0
      }
      if (!"d_parasitaemia" %in% names(self$population)) {
        self$population$d_parasitaemia = 0
      }
      
      self$locations = tibble::tibble(locations, gametocyte_load = NA_real_, EIR = NA_real_)
      self$mosquito = mosquito
      self$set_mosquito_raster(mosquito_raster)
      self$p_relapse = p_relapse
      self$mean_recovery = mean_recovery
      
      # Initialise log
      if ("linelist" %in% log_options) {
        self$log$linelist = DynamicFrame$new(
          ID = self$population$ID[!is.na(self$population$t_infection)],
          t_infection = self$population$t_infection[!is.na(self$population$t_infection)],
          source = "Seed"
        )
      }
      if ("compartment" %in% log_options) {
        self$log$compartment = DynamicFrame$new(
          t = self$t,
          infected = sum(!is.na(self$population$t_infection)))
      }
      if ("EIR" %in% log_options) {
        self$log$EIR = DynamicFrame$new(
          t = self$t,
          self$vis_grid,
          ento_inoculation_rate = apply(self$vis_grid, 1, self$calculate_EIR)
        )
      }
      
      invisible(self)
    },
    
    
    #' @description Iterate simulation
    #' @details
    #' Advance simulation by `dt`. Calculates the current state of population, infects
    #' mosquitoes and population simultaneously, updates population and mosquitoes, then
    #' stores observations for analysis.
    #' 
    #' @param dt Time to advance simulation by in days
    #' @param debug Optional, print debug information
    iterate = function(dt, debug=FALSE) {
      self$t = self$t + dt
      self$min_dt = min(self$min_dt, dt)
      
      ## Calculate current state of human infectivity and immunity
      # If a human is due to recover, wipe the infection
      recovery_IDs = which(self$population$t_recovery <= self$t)
      self$population$t_infection[recovery_IDs] = NA # wipe infection
      self$population$t_recovery[recovery_IDs] = Inf # unschedule recovery
      # Otherwise calculate their disease state
      self$population$p_blood_gametocyte = self$human$infectivity(self$population, self$t)
      self$population$immunity = self$human$immunity(self$population, self$t)
      
      # Sum total human gametocyte load per location
      self$locations$gametocyte_load = 0
      infected = self$population[self$population$p_blood_gametocyte > 0,]
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
        self$t - self$mosquito_infections$t_inoculation <= self$mosquito$max_lifespan,]
      # Calculate current infectivity of infected mosquito clouds integrated over space
      # i.e., number of live mosquitoes with mature sporozoites from this event
      self$mosquito_infections$density = self$mosquito_infections$count *
        self$mosquito$survival(self$t - self$mosquito_infections$t_inoculation) *
        self$mosquito$sporogony(self$t - self$mosquito_infections$t_inoculation)
      
      # Expose locations
      self$locations$EIR = apply(self$locations, 1, self$calculate_EIR)
      
      # Expose and infect population who are not immune or are due to relapse
      susceptible = self$population[self$population$immunity < 1 | self$population$t_relapse <= self$t,]
      susceptible_EIR = numeric() # otherwise it'll be undefined
      if (nrow(susceptible) > 0) {
        susceptible_EIR = apply(
          susceptible,
          1,
          function(h) {
            ix = h$location_ix
            # get weighted average of EIR at the person's locations
            sum(self$locations$EIR[ix] * h$location_proportions)
          }
        ) * (1 - susceptible$immunity)
      }
      # shortcut for P(x > 0), x ~Pois(EIR*dt) OR relapse
      infect_ix = runif(nrow(susceptible)) > exp(-susceptible_EIR * dt) |
        susceptible$t_relapse <= self$t
      is_relapse = (susceptible$t_relapse <= self$t)[infect_ix]
      infect_IDs = susceptible$ID[infect_ix]
      if(any(is.na(infect_IDs))) print(susceptible_EIR)
      
      ## Update human population
      self$population$t_infection[infect_IDs] = self$t
      # clear and reschedule relapses and recoveries
      self$population$t_relapse[infect_IDs] = self$t + self$schedule_relapse(length(infect_IDs))
      self$population$t_recovery[infect_IDs] = self$t + self$schedule_recovery(length(infect_IDs))
      
      ## Update mosquito population
      # Expose mosquitoes at locations that have both gametocytes and mosquitoes
      new_mosquito_infections = self$locations[self$locations$gametocyte_load > 0 &
                                                 self$locations$mosquito_density > 0,
                                               c("x", "y", "gametocyte_load", "mosquito_density")]
      if (nrow(new_mosquito_infections) > 0) {
        # Infect mosquitoes at rate (assume no depletion of the susceptible compartment)
        new_mosquito_infections$count = new_mosquito_infections$gametocyte_load *
          new_mosquito_infections$mosquito_density * # N.B. local density could potentially be recomputed
          self$mosquito$bite_rate *
          dt
        new_mosquito_infections$density = 0 # gets overwritten later anyway but prevent NAs
        new_mosquito_infections$t_inoculation = self$t
        new_mosquito_infections$gametocyte_load <- new_mosquito_infections$mosquito_density <- NULL
        # Add new potential mosquitoes to list
        self$mosquito_infections = bind_rows(self$mosquito_infections, new_mosquito_infections)
      }
      
      # Log observations
      if ("linelist" %in% names(self$log)) {
        self$log$linelist$append(
          ID = infect_IDs,
          t_infection = self$t,
          source = as.character(ifelse(is_relapse, "Relapse", "Transmission"))
        )
      }
      if ("compartment" %in% names(self$log)) {
        self$log$compartment$append(
          t = self$t,
          infected = sum(!is.na(self$population$t_infection))
        )
      }
      if ("EIR" %in% names(self$log)) {
        self$log$EIR$append(
          t = self$t,
          self$vis_grid,
          ento_inoculation_rate = apply(self$vis_grid, 1, self$calculate_EIR)
        )
      }
      
      invisible(self)
    },
    
    #' @description Recompute location mosquito densities from a raster
    #' @details Updates simulation's raster and densities at each location
    #' 
    #' @param mosquito_raster Raster of mosquito densities
    set_mosquito_raster = function(mosquito_raster) {
      self$locations$mosquito_density = my_extract(mosquito_raster, self$locations[c("x", "y")])
      self$mosquito_raster = mosquito_raster
      invisible(mosquito_raster)
    },
    
    #' @description Entomological inoculation rate
    #' @details Calculate EIR at a location with x and y elements
    #' 
    #' @param loc List or dataframe with `x` and `y` columns
    #' 
    #' @return Numeric vector of EIR
    calculate_EIR = function(loc) {
      # Calculate total EIR=HBR*SIR from all mosquito infection events
      with(
        self$mosquito_infections,
        sum(density *
              self$mosquito$bite_rate *
              # Disperse load over total cloud area (TODO: check integral=1)
              self$mosquito$migration(x-loc["x"], y-loc["y"], self$t-t_inoculation) *
              self$sporozoite_infection_rate)
      )
    },
    
    relapse_shape = 1,
    relapse_rate = 1/30,
    #' @description Human relapse distribution
    #' @details
    #' Schedule a relapse infection for number of days after the primary infection.
    #' Returns `Inf` value if no relapse.
    #' shape = (mean / sd)^2
    #' rate =  mean / sd^2
    #' 
    #' @param n Number of people to schedule relapses for
    schedule_relapse = function(n) {
      ifelse(runif(n) < self$p_relapse,
             14 + rgamma(n, self$relapse_shape, self$relapse_rate),
             Inf)
    },
    
    #' @description Human recovery distribution
    #' @details
    #' Schedule a recovery for number of days after the primary infection. Can
    #' also be used to schedule the success of treatment.
    #' Returns `Inf` value if no complete clearance
    #' 
    #' @param n Number of people to schedule relapses for
    schedule_recovery = function(n) {
      # rexp(n, rate = 1 / self$mean_recovery)
      Inf
    },
    
    #' @description Plot current state of a simulation
    #' @details
    #' Plots map with epicurve or the starting parameters if simulation is in its initialised state
    #' 
    #' @param x Simulation object
    #' @param t Optional, plot simulation as it was at this time
    #' @param init Logical, show initialisation state with mosquito raster
    #' @param ... Additional arguments
    plot = function(t=NULL, init=NULL, ...) {
      if (is.null(init)) {
        # check if sim is still in its initial state
        init = is.infinite(self$min_dt)
      }
      if (init) {
        plot_state(self, t, background="mosquito", ...)
      } else {
        p_state = plot_state(self, t, ...)
        p_epicurve = plot_epicurve(self, t, ...) + guides(fill = "none")
        p_state / p_epicurve + plot_layout(heights=c(5,1), guides="collect")
      }
    },
    
    #' @description Print simulation
    #' 
    #' @param x Simulation object
    #' @param ... Additional arguments
    print = function(...) {
      print(
        list(
          "t" = self$t,
          "linelist" = self$linelist,
          "population" = self$population,
          "attack_rate" = length(unique(self$linelist$ID)) / nrow(self$population)
        )
      )
    }
  ),
  
  
  #' @field population_collapse Represent population with only one set of coordinates
  #' @field population_expand Represent population at all of their coordinates
  #' @field linelist Safely access linelist if logged
  #' @field compartment Safely access compartments if logged
  #' @field EIR Entomological inoculation rate by coordinates over time
  #' @field vis_grid Data frame version of `vis_raster`
  #' @field attack_rate Proportion of population who have been infected at least once
  active = list(
    
    linelist = function() {
      stopifnot("Tried to access linelist logs but 'linelist' not specified in `log_options`" = "linelist" %in% names(self$log))
      self$log$linelist$df
    },
    
    compartment = function() {
      stopifnot("Tried to access compartment size logs but 'compartment' not specified in `log_options`" = "compartment" %in% names(self$log))
      self$log$compartment$df
    },
    
    EIR = function() {
      stopifnot("Tried to access EIR logs but 'EIR' not specified in `log_options`" = "EIR" %in% names(self$log))
      self$log$EIR$df
    },
    
    vis_grid = function() {
      as.data.frame(self$vis_raster, xy=T)[, c("x", "y")]
    },
    
    #' @description Representation of population with only one set of coordinates
    #' @details Just places people at their highest proportion location
    population_collapse = function() {
      self$population %>%
        mutate(location_ix = unlist(map(location_ix, first)),
               x = self$locations$x[location_ix],
               y = self$locations$y[location_ix])
    },
    
    #' @description Represent population at all of their coordinates
    population_expand = function() {
      self$population %>%
        unnest(cols=c(location_ix, location_proportions)) %>%
        mutate(x = self$locations$x[location_ix],
               y = self$locations$y[location_ix])
    },
    
    #' @description Fraction of non-seeded people who were infected
    attack_rate = function() {
      n_seeds = sum(self$linelist$source == "Seed")
      n_infections = length(unique(self$linelist$ID))
      n_people = nrow(self$population)
      (n_infections - n_seeds) / (n_people - n_seeds)
    }
  )
)
