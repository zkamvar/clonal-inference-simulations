#' Filter loci that don't have the correct minor allele frequencies
#'
#' @param x a genlight or snpclone object
#' @param maf the minor allele frequency (default to 0.05)
#' @param e an error term added to maf
#'
#' @return a genlight or snpclone object
#' @export
#'
#' @examples
#' ex <- structure(c(1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L,
#'       0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'       2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 0L, 0L, 0L, 0L, 0L,
#'       0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 2L,
#'       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L,
#'       2L, 1L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
#'       0L, 0L, 0L, 1L, 0L, 0L, 0L), .Dim = c(10L, 10L), .Dimnames = list(
#'           c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), c("0_0",
#'           "0_1", "0_2", "0_3", "0_4", "0_5", "0_6", "0_7", "0_8", "0_9"
#'           )))
#' gl <- new("genlight", ex)
#' gl
#' filter_maf(gl)
filter_maf <- function(x, maf = 0.05, e = .Machine$double.eps){
  x[, mafx(as.matrix(x), maf + e)]
}

#' Indicate if loci have the correct MAF
#'
#' @param x a matrix derived from a genlight object
#' @param maf minor allele frequency
#'
#' @return a vector of logical values
#' @noRd
#'
#' @examples
#' ex <- structure(c(1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 0L, 1L, 0L,
#'       0L, 0L, 0L, 0L, 0L, 0L, 0L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
#'       2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 0L, 0L, 0L, 0L, 0L,
#'       0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 2L,
#'       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 2L, 1L, 2L,
#'       2L, 1L, 2L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 0L, 0L,
#'       0L, 0L, 0L, 1L, 0L, 0L, 0L), .Dim = c(10L, 10L), .Dimnames = list(
#'           c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), c("0_0",
#'           "0_1", "0_2", "0_3", "0_4", "0_5", "0_6", "0_7", "0_8", "0_9"
#'           )))
#' mafx(ex, 0.051)
mafx <- function(x, maf){
  s <- colSums(x, na.rm = TRUE)/(2*nrow(x));
  return(s >= maf & (1 - s) >= maf)
}
