#' Generate a tidy table of diversity statistics
#'
#' @param gid a genind or genclone object
#' @param n number of bootstrap replicates (default: 1000)
#' @param ... arguments to be passed on to \code{\link[poppr]{diversity_stats}}
#' @param verbose when \code{TRUE}, a status message will be printed at the
#'   beginning of the analysis. \code{FALSE} suppresses this message.
#' @return a tidy data frame with estimates and variances
#' @export
#'
#' @examples
#' uSimp <- function(x){
#'   lambda <- vegan::diversity(x, "simpson")
#'   x <- drop(as.matrix(x))
#'   if (length(dim(x)) > 1){
#'     N <- rowSums(x)
#'   } else {
#'     N <- sum(x)
#'   }
#'   return((N/(N-1))*lambda)
#' }
#' library("poweRlaw")
#' power_law_beta <- function(x){
#'   if (length(x) == 1) return(NA_real_)
#'   xpow <- displ(x[x > 0])                 # Generate the distribution
#'   xpow$setPars(estimate_pars(xpow))       # Estimate the parameters
#'   xdat <- plot(xpow, draw = FALSE)        # Extract the data
#'   xlm  <- lm(log(y) ~ log(x), data = xdat) # Run log-log linear model for slope
#'   return(-coef(xlm)[2])
#' }
#'
#' Beta <- function(x){
#'   x <- drop(as.matrix(x))
#'   if (length(dim(x)) > 1){
#'     res <- apply(x, 1, power_law_beta)
#'   } else {
#'     res <- power_law_beta(x)
#'   }
#'   return(res)
#' }
#' library("poppr")
#' data(Pinf)
#' seppop(Pinf) %>%
#'   lapply(tidy_diversity, n = 100, B = Beta, uSimp = uSimp) %>%
#'   bind_rows()
tidy_diversity <- function(gid, n = 1000, ..., verbose = TRUE){
  tab       <- mlg.table(gid, plot = FALSE)
  the_pops  <- popNames(gid)
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
  the_stats <- poppr::diversity_stats(tab, ...)
  if (!is.matrix(the_stats)){
    cnames              <- names(the_stats)
    dim(the_stats)      <- c(1, length(the_stats))
    colnames(the_stats) <- cnames
  }
  the_stats <- the_stats %>%
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
