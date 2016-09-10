#' Apply a function on a single element list
#'
#' This allows one to work with list columns in dplyr
#'
#' @param x a single element list
#' @param f a function to apply to the first element of x
#'
#' @return the output of f
#' @export
#'
#' @examples
#' x <- list(y = 1:10)
#' listfun(x, sum)
listfun <- function(x, f){
  f <- match.fun(f)
  f(x[[1]])
}
