#' Gather statistics on locus contribution
#'
#' @param gid a genclone object
#' @param n the number of permutations for calculating the genotype accumulation curve.
#' @param verbose should this be noisy?
#' @param ... arguments to be passed on to diversity_stats.
#'
#' @return a tibble object
#' @export
#'
#' @examples
#' data(monpop)
#' mll(monpop) <- "original"
#' monpop <- monpop %>%
#'   splitStrata(~Tree/Year/Symptom) %>%
#'   setPop(~Tree)
#' monpop %>%
#'   seppop() %>%
#'   lapply(tidy_locus_contribution) %>%
#'   bind_rows()
tidy_locus_contribution <- function(gid, n = 999, verbose = TRUE, ...){
  if (verbose){
    msg <- "Calculating locus diversity for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  NMLL     <- nmll(gid)
  rrmat    <- rrmlg(gid)
  pgen_gid <- poppr::pgen(gid, mul = 1/2, log = FALSE, by_pop = FALSE)
  gc       <- genotype_curve(gid, plot = FALSE, sample = n, quiet = TRUE)
  cont_div <- locus_contribution(rrmat, NMLL) %>%
    diversity_stats(...) %>%
    as.list()
  out <- c(cont_div,
           list(rrmat = ~list(rrmat),
                curve = ~list(gc),
                pgen  = ~list(pgen_gid),
                pop   = ~popNames(gid))
           )
  return(data_frame_(out))
}


#' Compute the number of unique multilocus genotypes contributed by each locus.
#'
#' @param rrmat a matrix of round-robin multilocus genotypes for each locus
#' @param NMLL the number of multilocus genotypes in the sample
#'
#' @return
#' @noRd
#'
#' @examples
#' data(monpop)
#' mll(monpop) <- "original"
#' locus_contribution(rrmlg(monpop), nmll(monpop))
locus_contribution <- function(rrmat, NMLL){
  NMLL - colSums(!apply(rrmat, 2, duplicated), na.rm = TRUE)
}
