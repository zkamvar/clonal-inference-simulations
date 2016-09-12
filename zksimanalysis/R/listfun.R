#' Apply a function on a single element list
#'
#' This allows one to work with list columns in dplyr
#'
#' @param x a single element list
#' @param f a function to apply to the first element of x
#' @param ... any arguments that go with f
#'
#' @return the output of f
#' @export
#'
#' @examples
#' x <- list(y = c(1:10, NA))
#' listfun(x, sum, na.rm = TRUE)
listfun <- function(x, f, ...){
  f <- match.fun(f)
  f(x[[1]], ...)
}
