#' Generate a tidy table of diversity statistics
#'
#' @param gid a genind or genclone object
#' @param n number of bootstrap replicates (default: 1000)
#' @param ... arguments to be passed on to \code{\link[poppr]{diversity_stats}}
#' @param verbose
#' @param strata
#'
#' @return
#' @export
#'
#' @examples
tidy_diversity <- function(gid, n = 1000, ..., verbose = TRUE, strata = NA){
  tab       <- mlg.table(gid, plot = FALSE)
  the_pops  <-  popNames(gid)
  the_pops  <- if (is.null(the_pops)) NA_character_ else the_pops
  if (verbose){
    msg <- "Calculating diversity for "
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  the_pops  <- data_frame(pop = the_pops)
  the_stats <- poppr::diversity_stats(tab, ...) %>%
    as_data_frame() %>%
    bind_cols(the_pops, .)
  the_boot  <- diversity_boot(tab, n = n, ...)
  the_vars  <- get_boot_var(the_boot) %>%
    as_data_frame()
  return(bind_cols(the_stats, the_vars))
}

#' Function modified from poppr's get_boot_se
#'
#' @param bootlist a list of boot objects
#'
#' @return a matrix
#' @noRd
#'
get_boot_var <- function(bootlist){
  npop     <- length(bootlist)
  bstats   <- bootlist[[1]]$t0
  nstat    <- length(bstats)
  matnames <- list(Pop = names(bootlist), Index = names(bstats))
  THE_FUN  <- function(x) apply(x$t, 2, stats::var, na.rm = TRUE)
  resmat   <- matrix(nrow = npop, ncol = nstat, dimnames = matnames)
  resmat[] <- t(vapply(bootlist, FUN = THE_FUN, FUN.VALUE = bstats))
  colnames(resmat) <- paste0(colnames(resmat), ".var")
  return(resmat)
}
