library(R6)

#' Dynamic expanding data frame class
#' 
#' Designed to double space in memory when needed instead of always expanding the list.
#' Dynamically assembles the resulting data frame from a list of elements when called.
#' 
#' This is good if we:
#' 1. Don't know the total data size ahead of time
#' 2. Add rows incrementally and often, e.g., in a loop
#' 3. Retrieve the full data infrequently, e.g., while plotting.
DynamicFrame = R6Class(
  "DynamicFrame",
  
  public = list(
    #' @description Initialise a data frame with some data
    #' @param ... Add arguments as rows of a dataframe
    initialize = function(...) {
      self$append(...)
    },
    
    #' @description Add rows to the data
    #' @param ... Data as if passed to `tibble::tibble(...)`
    append = function(...) {
      if (private$size == private$count) {
        length(private$data) = private$size = private$size * 2
      }
      private$count = private$count + 1
      private$data[[private$count]] = tibble::tibble(...)
      
      invisible(self)
    }
  ),
  
  #' @field data List of data
  #' @field size Current number of slots in memory
  #' @field count Number of elements stored
  private = list(
    data = list(),
    size = 0,
    count = 0
  ),
  
  #' @field df Retrieve data as a tibble
  active = list(
    df = function() {
      tibble::as_tibble(data.table::rbindlist(private$data))
    }
  )
)
