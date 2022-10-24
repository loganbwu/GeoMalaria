library(R6)

#' Shell class to hold human attributes and methods
#' 
#' To use a different mosquito, either:
#' 1. Use `$new(...)$` to initialise an instance of this class with different constants
#' 2. Use the `$set(...)` method to overwrite class methods, then initialise an instance
#' 3. Create a child class
Human = R6Class(
  "Human",
  
  #' @field duration_infectivity Constrain infectious period
  public = list(
    duration_infectivity = NULL,
    
    #' @description Initialise a mosquito class, optionally overwriting the default parameters
    #' 
    #' @param ... Named arguments of parameter constants or functions
    initialize = function(duration_infectivity = 100) {
      self$duration_infectivity = duration_infectivity
    },
    
    #' @description Lifecycle stages for P. vivax in population after inoculation
    #' @details
    #' Calculates the blood parasite density after infection in counts per microlitre.
    #' Specify either `population` AND `t`, OR `dt`
    #' 
    #' @param population Dataframe of population with a `t_infection` column
    #' @param t Simulation time
    #' @param dt Time since infection
    parasitaemia = function(population=NULL, t=NULL, dt=NULL) {
      stopifnot("Function needs (population AND t) OR dt" = xor(!is.null(population) & !is.null(t), !is.null(dt)))
      if (is.null(dt)) {
        dt = t - population$t_infection
      }
      
      parasite_density = approx(c(12, 20, 40, self$duration_infectivity),
                                c(0, 1e9, 1e9, 0),
                                dt,
                                yleft = 0, yright = 0)$y
      replace_na(parasite_density, 0)
    },
    
    #' @description Thresholds for human clinical outcomes
    #' 
    #' @param d_parasitaemia
    #' @name human_clinical
    
    #' @description Criteria to show symptoms
    #' 
    #' @rdname human_clinical
    has_symptoms = function(d_parasitaemia) {
      d_parasitaemia > 1e8
    },
    
    #' @description Criteria for death
    #'
    #' @rdname human_clinical
    has_death = function(d_parasitaemia) {
      d_parasitaemia > 5e8
    },
    
    #' @description Lifecycle stages for P. vivax in population after inoculation
    #' @details
    #' Calculates a vector of probability of infecting mosquito with gametocytes in the grid cell.
    #' Specify either `population` AND `t`, OR `dt`
    #' 
    #' @param population Dataframe of population with a `t_infection` column
    #' @param t Simulation time
    #' @param dt Time since infection
    infectivity = function(population=NULL, t=NULL, dt=NULL) {
      stopifnot("Function needs (population AND t) OR dt" = xor(!is.null(population) & !is.null(t), !is.null(dt)))
      if (is.null(dt)) {
        dt = t - population$t_infection
      }
      
      infectivity = approx(c(12, 17, 20, self$duration_infectivity),
                           c(0, 1, 1, 0),
                           dt,
                           yleft = 0, yright = 0)$y
      replace_na(infectivity, 0)
    },
    
    #' @description Human immunity profile over time
    #' @details
    #' Full immunity then decays to 50% 
    #' scales # bites not # susceptible people. E.g. no one retains complete immunity
    #' after decay, everyone retains partial protection.
    #' 
    #' Specify either `population` AND `t`, OR `dt`
    #' 
    #' @param population Dataframe of population with a `t_infection` column
    #' @param t Simulation time
    #' @param dt Time since infection
    immunity = function(population=NULL, t=NULL, dt=NULL) {
      stopifnot("Function needs (population AND t) OR dt" = xor(!is.null(population) & !is.null(t), !is.null(dt)))
      if (is.null(dt)) {
        dt = t - population$t_infection
      }
      
      immunity = approx(c(0, self$duration_infectivity, 2*self$duration_infectivity),
                        c(1, 1, 0),
                        dt,
                        yleft = 0, yright = 1)$y
      replace_na(immunity, 0)
    },
    

    #' @description Plot characteristics
    #' @examples
    #' human = Human$new()
    #' plot(human)
    plot = function() {
      plot_human_infectivity(self) / plot_human_clinical_outcomes(self)
    }
  )
)


#' Plot human infectivity curve
#' 
#' @param human Human object
plot_human_infectivity = function(human) {
  # Add visible bindings
  value <- name <- NULL
  
  data = tibble(dt = seq(0, 2*human$duration_infectivity, length.out=1000),
         gametocyte_load = human$infectivity(dt=dt),
         immunity = human$immunity(dt=dt))
  
  data %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name)) +
    geom_line() +
    scale_y_continuous(labels = label_auto) +
    labs(title = "Inoculation load and immunity after infection",
         x = "Days since human infection", y = "Probability", color = "Attribute")
}

plot_human_clinical_outcomes = function(human) {
  
  data = tibble(dt = seq(0, 2*human$duration_infectivity, length.out=1000),
                parasitaemia = human$parasitaemia(dt=dt))
  cutoffs = tribble(
    ~t, ~name,
    min(data$dt[human$has_symptoms(data$parasitaemia)]), "Symptoms",
    min(data$dt[human$has_death(data$parasitaemia)]), "Death"
  )
  
  ggplot(data, aes(x = dt, y = parasitaemia)) +
    geom_line() +
    geom_vline(data=cutoffs, aes(xintercept=t, color=name), show.legend=F) +
    geom_text(data=cutoffs, aes(x=t, y=0, label=name, color=name), angle=90, hjust=0, vjust=1.2, show.legend=F) +
    scale_y_continuous(labels = label_auto) +
    labs(title = "Clinical timeline",
         x = "Days since human infection", y = "Parasitaemia")
}

# hum = Human$new()
# plot(hum)
