#' Generate a tidy data frame of the index of assocation
#'
#' @param gid a genind or genclone object
#' @param ... any parameters to be passed to \code{\link[poppr]{ia}}
#' @param verbose when \code{TRUE}, a status message will be printed at the
#'   beginning of the analysis. \code{FALSE} suppresses this message.
#' @param keepdata when \code{TRUE}, the data set is kept in a column called
#'   "dataset"
#'
#' @return a tibble from the ia object
#' @export
#'
#' @examples
#' library('poppr')
#' data(nancycats)
#' nancycats %>%
#'   seppop() %>%
#'   lapply(tidy_ia, sample = 99, hist = FALSE) %>%
#'   bind_rows
tidy_ia <- function(gid, ..., verbose = TRUE, keepdata = TRUE){
  if (verbose){
    msg <- "Calculating ia for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  res  <- poppr::ia(gid, ..., valuereturn = TRUE)
  pops <- popNames(gid)
  pops <- if (length(pops) > 1) list(pops) else pops
  if (inherits(res, "ialist")){
    vals <- res[[1]]
    out  <- list(Ia         = ~vals["Ia"],
                 p.Ia       = ~vals["p.Ia"],
                 rbarD      = ~vals["rbarD"],
                 p.rD       = ~vals["p.rD"],
                 samples.ia = ~list(res[[2]]$Ia),
                 samples.rd = ~list(res[[2]]$rbarD),
                 pop        = ~pops
                 )
  } else {
    out <- list(Ia    = ~res["Ia"],
                rbarD = ~res["rbarD"],
                pop   = ~pops
                )
  }
  if (keepdata){
    out <- c(out, list(dataset = ~list(gid)))
  }
  return(data_frame_(out))
}
