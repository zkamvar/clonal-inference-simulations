#' Give a tidy locus table
#'
#' @param gid a genclone object
#' @param ... parameters to be passed on
#' @param verbose shoudl this be noisy?
#'
#' @return a tibble
#' @export
#'
#' @examples
#' data(nancycats)
#'   seppop(nancycats) %>%
#'   lapply(tidy_locus_table) %>%
#'   bind_rows()
tidy_locus_table <- function(gid, ..., verbose = TRUE){
  if (verbose){
    msg <- "Calculating locus diversity for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  tab <- locus_table(gid, information = FALSE, ...)
  avg <- tab["mean", ] %>% as.list()
  out <- list(tab = ~list(tab[-nrow(tab), ]), pop = ~popNames(gid))
  out <- c(avg, out)
  return(data_frame_(out))
}
