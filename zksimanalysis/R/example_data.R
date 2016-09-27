#' A wrapper to get example data
#'
#' @return a string with the path to the example
#' @export
#'
#' @examples
#' example_data()
example_data <- function(){
  system.file("files", "sample_data.feather", package = "zksimanalysis")
}

#' Example simulation data
#'
#' @name sample_data
#' @docType data
#' @usage data(sample_data)
#' @description This data set contains simulations of clonal populations
#' @format a tibble containing genclone objects
NULL
