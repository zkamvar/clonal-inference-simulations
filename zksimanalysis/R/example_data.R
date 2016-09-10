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
