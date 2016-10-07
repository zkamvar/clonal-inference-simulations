#' Give a tidy locus table
#'
#' @param gid a genclone object
#' @param ... parameters to be passed on
#' @param verbose should this be noisy?
#' @param cc should the data be clone corrected as well?
#' @param strata The population strata to use with cc.
#'
#' @return a tibble
#' @export
#'
#' @examples
#' data(nancycats)
#'   seppop(nancycats) %>%
#'   lapply(tidy_locus_table) %>%
#'   bind_rows()
tidy_locus_table <- function(gid, ..., verbose = TRUE, cc = TRUE, strata = NA){
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
  out <- list(tab = ~list(tab[-nrow(tab), ]))
  if (cc){
    tabcc <- gid %>% 
      clonecorrect(strata = strata) %>% 
      locus_table(information = FALSE, ...)
    avgcc        <- tabcc["mean", ] %>% as.list()
    names(avgcc) <- paste0(names(avgcc), "cc")
    out          <- c(out, avgcc, list(tabcc = ~list(tabcc[-nrow(tabcc), ])))
  }
  out <- c(avg, out, list(pop = ~popNames(gid)))
  return(data_frame_(out))
}
