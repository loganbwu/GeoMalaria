library(R6)

#' Relapse scheduler for P. vivax
Relapse = R6Class(
  "Relapse",
  public = list(
    
    schedule = tibble(ID = integer(), t = numeric()),
    
    initialize = function() {
      
    },
    
    print = function() {
      print(self$schedule)
    },
    
    schedule_stochastic_relapse = function(ID) {
      self$schedule_relapse(ID, runif(length(ID), 60, 80))
      invisible(self$schedule)
    }
    
    schedule_relapse = function(ID, t) {
      self$schedule = bind_rows(self$schedule, tibble(ID, t))
      invisible(self$schedule)
    },
    
    #' Get IDs of people due to relapse and remove from the schedule
    get_relapses = function(t) {
      relapses = self$schedule$ID[self$schedule$t <= t]
      self$schedule = schedule[self$schedule$t > t,]
      return(relapses)
    }
  )
)