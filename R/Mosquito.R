library(R6)

#' Shell class to hold species attributes and methods
#' 
#' To use a different mosquito, either:
#' 1. Use `$new(...)$` to initialise an instance of this class with different constants
#' 2. Use the `$set(...)` method to overwrite class methods, then initialise an instance
#' 3. Create a child class
Mosquito = R6Class(
  "Mosquito",
  
  #' @field distribution Raster of spatial mosquito counts
  #' @field death_rate Proportion of mosquitoes that die per day, simulated continuously
  #' @field travel Constant of mosquito travel in variance of location after one day
  #' @field bite_rate Bites per mosquito per density
  #' @field max_lifespan Automatic: ignore mosquitoes after this duration after biting
  #' @field max_flight_range Automatic: ignore mosquitoes this distance from original bite site
  public = list(
    distribution = NULL,
    death_rate = NULL,
    travel = NULL,
    bite_rate = NULL,
    max_lifespan = 0,
    max_flight_range = NULL,
    
    #' @description Initialise a mosquito class, optionally overwriting the default parameters
    #' 
    #' @param ... Named arguments of parameter constants or functions
    initialize = function(distribution,
                          death_rate = 0.25,
                          travel = 1^2,
                          bite_rate = 0.1) {
      self$death_rate = death_rate
      self$travel = travel
      self$bite_rate = bite_rate
      
      # Derive bound for the vast majority of the mosquito lifespan
      while (self$survival(self$max_lifespan) > 0.01) {
        self$max_lifespan = self$max_lifespan + 1
      }
      
      # Derive max range to consider mosquito diffusion
      self$max_flight_range = qnorm(0.975, sd=sqrt(self$travel * self$travel))
      
    },
    
    #' @description Proportion of mosquitoes that survive over time
    #' 
    #' @param dt Time since mosquito inoculation in days
    survival = function(dt) {
      survival = dexp(dt, rate=self$death_rate) / self$death_rate
    },
    
    #' @description Proportion of mosquitoes with sporozoites over time
    #' 
    #' @param dt Time since mosquito inoculation in days
    sporogony = function(dt) {
      approx(c(8, 10), c(0, 1), dt, yleft=0, yright=1)$y
    },
    
    #' @description Mosquito density
    #' @details
    #' Mosquito density calculated if they started at a point and traveled `dx` and
    #' `dy` in time `dt`. Variance is constant at (km/day)^2. Uses the Laplace isometric
    #' diffusion solution.
    #' 
    #' @param dx Distance traveled horizontally in km
    #' @param dy Distance traveled vertically
    #' @param dt Time since origin in days
    migration = function(dx, dy, dt) {
      exp(dnorm(dx, sd=sqrt(self$travel * dt), log=T) + 
            dnorm(dy, sd=sqrt(self$travel * dt), log=T))
    },
    
    #' @description Print mosquito summary
    #' 
    
    #' @description Plot mosquito characteristics
    #' @examples
    #' mos = Mosquito$new()
    #' plot(mos)
    plot = function() {
      plot_mosquito_migration(self) / plot_mosquito_infectivity(self)
    }
  )
)


#' Plot spread of mosquitoes over distance
#' 
#' @param mos Mosquito object
#' @examples
#' mos = Mosquito$new()
#' plot_mosquito_migration(mos)
plot_mosquito_migration = function(mos) {
  dx = seq(-mos$max_flight_range,
           mos$max_flight_range,
           length.out = 100)
  dt = seq(1, mos$max_lifespan, length.out=7)
  expand.grid(dx=dx, dt=dt) %>%
    mutate(density = mos$migration(dx, 0, dt)) %>%
    ggplot(aes(x=dx, y=density, color=dt, group=dt)) +
    geom_line() +
    scale_y_continuous(labels = label_auto) +
    scale_color_binned(breaks = dt) +
    labs(title = "Distance travelled by mosquitoes",
         x = "Distance (1D shown; actual diffusion is 2D)", y = "Density",
         color = "Days since\nblood meal")
}

#' Plot mosquito infectivity curve
#' 
#' @param mos Mosquito object
#' @examples
#' mos = Mosquito$new()
#' plot_mosquito_infectivity(mos)
plot_mosquito_infectivity = function(mos) {
  # Add visible bindings
  survival <- sporozoites <- value <- name <- NULL
  
  tibble(dt = seq(0, mos$max_lifespan, length.out=1000),
         survival = mos$survival(dt),
         sporozoites = mos$sporogony(dt),
         infectivity = survival*sporozoites) %>%
    pivot_longer(cols = -dt) %>%
    ggplot(aes(x = dt, y = value, color = name, linetype = name)) +
    geom_line() +
    scale_y_continuous(labels = label_auto) +
    labs(title = "Inoculation probability after blood meal",
         x = "Days since blood meal", y = "Probability", color = "Attribute", linetype = "Attribute")
}
