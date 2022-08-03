library(R6)

Simulation = R6Class(
  "Simulation",
  public = list(
    # Constants
    t = 0,
    mosquito_raster = NULL,
    n_people = NULL,
    max_mosquito_lifespan = NULL,
    max_human_infectivity = NULL,
    max_mosquito_flight_range = NULL,
    infect_prob_unit = NULL,
    
    # States
    people = NULL,
    mosquito_infections = tibble(),
    people_infections = tibble(),
    
    initialize = function(env_dimensions,
                          n_people,
                          max_mosquito_lifespan, # days
                          max_human_infectivity, # days
                          max_mosquito_flight_range, # km
                          infect_prob_unit) { # probability of infection for 1 mosquito/cell/day
      
      # Raster for environmental effects
      private$env_raster = raster(vals=0, nrows=100, ncols=100,
                                  xmn=0, xmx=env_dimensions[1], ymn=0, ymx=env_dimensions[2], crs=NA)
      # Raster for fine visualisation of continuous fields
      private$vis_raster = raster(vals=0, nrows=100, ncols=100,
                                  xmn=0, xmx=env_dimensions[1], ymn=0, ymx=env_dimensions[2], crs=NA)
      
      # Parameters
      self$n_people = n_people
      self$max_mosquito_lifespan = max_mosquito_lifespan
      self$max_human_infectivity = max_human_infectivity
      self$max_mosquito_flight_range = max_mosquito_flight_range
      self$infect_prob_unit = infect_prob_unit
      
      # Initialise data
      # People
      self$people = tibble(
        ID = seq_len(n_people),
        X = runif(n_people, max=env_dimensions[1]),
        Y = runif(n_people, max=env_dimensions[2]),
        infected = runif(n_people) < 0.1,
        t_infection = ifelse(infected, self$t, NA),
        blood_gametocyte_prob = 0
      )
      # Mosquitos
      self$mosquito_raster = private$env_raster
      noise = noise_perlin(dim(private$env_raster)[1:2])
      noise = (noise - min(noise)) * 1000
      self$mosquito_raster[] = noise
      names(self$mosquito_raster) = "count"
    },
    
    plot_init = function() {
      mosquito_data = as.data.frame(self$mosquito_raster, xy=T) %>%
        rename(X = x, Y = y)
      humans = self$people %>%
        mutate(Infection = case_when(t_infection == self$t ~ "Today",
                                     t_infection < self$t ~ "Historical",
                                     TRUE ~ "None"))
      ggplot(mapping = aes(x=X, y=Y)) +
        geom_raster(data = mosquito_data, aes(fill=count)) +
        geom_point(data = humans, aes(color=Infection)) +
        coord_equal() +
        scale_fill_viridis_c(option = "magma") +
        scale_color_manual(values = c("Today"="tomato",
                                      "Historical"="steelblue",
                                      "None"="white"))
    },
    
    iterate = function(dt) {
      self$t = self$t + dt
      # convert daily unit probability to probability across dt (units of days)
      k = -log(1-self$infect_prob_unit)
      log_S_dt = -k * dt
      
      # Calculate current state of human infectivity
      self$people$blood_gametocyte_prob = private$human_infectivity(self$t - self$people$t_infection)
      
      # Expose mosquitos and add to current infected mosquitos
      new_mosquito_infections = self$people %>%
        filter(blood_gametocyte_prob > 0) %>%
        # Infect mosquitos at rate (assume no susc. depletion)
        mutate(infected_biting_rate = blood_gametocyte_prob *
                 my_extract(self$mosquito_raster, .[c("X", "Y")]) *
                 dt,
               t_inoculation = self$t) %>%
        filter(infected_biting_rate > 0) %>%
        select(X, Y, infected_biting_rate, t_inoculation)
      
      # Add new events to existing events
      self$mosquito_infections = bind_rows(self$mosquito_infections,
                                           new_mosquito_infections) %>%
        # Remove expired events
        filter(self$t - t_inoculation < self$max_mosquito_lifespan) %>%
        # Calculate current infectivity of infected mosquito clouds integrated over space
        mutate(sporozoite_load = infected_biting_rate *
                 private$mosquito_infectivity(self$t - t_inoculation))
      
      # Expose and infect people
      new_human_infections = self$people %>%
        filter(!infected) %>%
        merge(self$mosquito_infections, by=character()) %>%
        as_tibble() %>%
        mutate(distance = sqrt((X.x-X.y)^2 + (Y.x-Y.y)^2)) %>%
        filter(distance < self$max_mosquito_flight_range) %>%
        mutate(ento_inoculation_rate = sporozoite_load *
                 # Disperse load over total cloud area (TODO: check integral=1)
                 private$mosquito_migration(distance, self$t-t_inoculation),
               log_noinfect_prob = ento_inoculation_rate * log_S_dt,
        ) %>%
        group_by(ID) %>%
        summarise(infect_prob = 1 - exp(sum(log_noinfect_prob))) %>%
        filter(runif(n()) < infect_prob)
      
      print(paste("Infected", nrow(new_human_infections), "people"))
      
      # Update population
      self$people$infected[new_human_infections$ID] = TRUE
      self$people$t_infection[new_human_infections$ID] = self$t
    },
    
    plot = function() {
      neg_k = log(1-self$infect_prob_unit)
      humans = self$people %>%
        mutate(Infection = case_when(t_infection == self$t ~ "New",
                                     t_infection < self$t ~ "Historical",
                                     TRUE ~ "None"))
      
      mosquitoes = as.data.frame(private$vis_raster, xy=T) %>%
        mutate(ID = row_number()) %>%
        rename(X = x, Y = y) %>%
        merge(self$mosquito_infections, by=character()) %>%
        as_tibble() %>%
        mutate(distance = sqrt((X.x-X.y)^2 + (Y.x-Y.y)^2)) %>%
        rename(X=X.x, Y=Y.x) %>%
        filter(distance < self$max_mosquito_flight_range) %>%
        mutate(ento_inoculation_rate = sporozoite_load *
                 private$mosquito_migration(distance, self$t-t_inoculation),
               log_noinfect_prob = ento_inoculation_rate * neg_k,
        ) %>%
        group_by(X, Y) %>%
        summarise(infect_prob = 1 - exp(sum(log_noinfect_prob)),
                  .groups = "drop")
      
      ggplot(mosquitoes, aes(x=X, y=Y)) +
        geom_raster(aes(fill=infect_prob), interpolate=TRUE) +
        scale_fill_viridis_c(option="inferno", limits=c(0, NA)) +
        new_scale("fill") +
        geom_point(aes(fill=blood_gametocyte_prob, color=Infection), data=humans, pch=21, size=3, stroke=2) +
        scale_fill_viridis_c(limits=c(0, 1)) +
        scale_color_manual(values=c("New"="tomato", "Historical"="steelblue", "None"=NA)) +
        coord_equal() +
        labs(title = paste("t =", self$t)) +
        theme(legend.key.height = unit(10, "pt"))
    }
  ),
  
  
  
  private = list(
    env_raster = NULL,
    vis_raster = NULL,
    
    #' Lifecycle stages for P. vivax in mosquitos after a blood meal
    #'
    #' @param t days since blood meal
    #' 
    #' @return boolean vector TRUE if mosquito is infectious
    vivax_mosquito_lifecycle = function(t) {
      t_sporogony = runif(length(t), 8, 10)
      infectious = t > t_sporogony
      ifelse(is.na(infectious), F, infectious)
    },
    
    #' Proportion of mosquitoes that survive over time
    mosquito_survival = function(dt) {
      rate = 0.5
      survival = dexp(dt, rate=rate) / rate
    },
    
    #' Infectivity of the vivax parasite in mosquitoes over time
    mosquito_sporogony = function(dt) {
      t_sporogony = runif(length(dt), 8, 10)
      as.numeric(dt > t_sporogony)
    },
    
    #' Overall infectivity of a mosquito after an infected blood meal
    mosquito_infectivity = function(dt) {
      private$mosquito_survival(dt) * private$mosquito_sporogony(dt)
    },
    
    # 1-D variance of mosquito travel in (km/day)^2
    mosquito_travel = 0.5^2,
    mosquito_migration = function(dx, dt) {
      dnorm(dx, sd=sqrt(private$mosquito_travel * dt))
    },
    
    #' Lifecycle stages for P. vivax in humans after inoculation
    #'
    #' @param t days since inoculation
    #' 
    #' @return vector of probability of infecting mosquito with gametocytes in the grid cell
    human_infectivity = function(dt) {
      t_incubation = runif(length(dt), 12, 17)
      infectivity = as.numeric(dt > t_incubation & dt < self$max_human_infectivity)
      replace_na(infectivity, 0)
    }
  )
)
