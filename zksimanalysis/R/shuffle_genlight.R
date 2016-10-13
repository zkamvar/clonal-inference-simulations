#' resample a genlight/genclone object and produce a
#'
#' @param x a snpclone or genlight object
#'
#' @return a snpclone or genlight object
#' @export
#'
#' @examples
#' gl <- glSim(100, 100, 100, parallel = FALSE)
#' gl
#' bitwise.ia(gl)
#' gls <- shuffle_genlight(gl)
#' gls
#' bitwise.ia(gls)
shuffle_genlight <- function(x){
  mat   <- as.matrix(x)
  shuff <- new("genlight",
               apply(mat, 2, sample),
               parallel = FALSE)
  x@gen <- shuff@gen
  return(x)
}

#' Calculate index of association from reshuffled genlight object
#'
#' @param i an integer, only used as a placeholder for the apply functions
#' @param x a genlight or snpclone object
#'
#' @return a value of the index of association
#' @export
#'
#' @examples
#' gl <- glSim(100, 100, 100, parallel = FALSE)
#' bitwise.ia(gl)
#' sample_bitwise_ia(1, gl)
sample_bitwise_ia <- function(i, x){
  bitwise.ia(shuffle_genlight(x))
}

#' A test for significance with genomic ia
#'
#' @param x a genlight or snpclone object
#' @param sample the number of samples for significance testing
#'
#' @return a list with the observed value and pvalue and the samples used.
#' @export
#'
#' @note This is a slow procedure. When tested on 185 samples over 10,000 SNPs,
#' this took about 5 seconds/sample, which is 84 minutes of computing time on
#' my mac.
#'
#' @examples
#' gl <- glSim(100, 200, 0, parallel = FALSE)
#' genomic_ia(gl, sample = 999)
genomic_ia <- function(x, sample = 0){
  obs  <- bitwise.ia(x)
  if (sample == 0) return(c(rbarD = obs))
  exp  <- vector(mode = "numeric", length = )
  exp  <- vapply(seq(sample), sample_bitwise_ia, numeric(1), x)
  pval <- (sum(exp >= obs) + 1)/(sample + 1)
  return(list(observed = c(rbarD = obs, p.rD = pval), samples = exp))
}
