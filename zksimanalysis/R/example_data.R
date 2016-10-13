#' A wrapper to get example data
#'
#' @param set the name of the data set you want to use. Defaults to "sample_data.feather"
#' @return a string with the path to the example
#' @export
#'
#' @examples
#' example_data()
example_data <- function(set = "sample_data.feather"){
  system.file("files", set, package = "zksimanalysis")
}

#' Example simulation data
#'
#' @name sample_data
#' @docType data
#' @usage data(sample_data)
#' @description This data set contains simulations of clonal populations
#' @format a tibble containing genclone objects
NULL
