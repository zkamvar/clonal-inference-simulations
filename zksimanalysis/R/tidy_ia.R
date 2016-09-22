#' Generate a tidy data frame of the index of association
#'
#' @param gid a genind or genclone object
#' @param ... any parameters to be passed to \code{\link[poppr]{ia}}
#' @param verbose when \code{TRUE}, a status message will be printed at the
#'   beginning of the analysis. \code{FALSE} suppresses this message.
#' @param keepdata when \code{TRUE}, the data set is kept in a column called
#'   "dataset"
#' @param cc return values of the clone-corrected data set (no re-sampling).
#' @param strata the strata at which to clone-correct. Defaults to entire data
#'   set.
#'
#' @return a tibble from the ia object
#' @export
#'
#' @examples
#' library('poppr')
#' data(partial_clone)
#' strata(partial_clone) <- data.frame(pop = pop(partial_clone))
#' partial_clone %>%
#'   seppop() %>%
#'   lapply(tidy_ia, sample = 999, plot = FALSE) %>%
#'   bind_rows()
tidy_ia <- function(gid, ..., verbose = TRUE, keepdata = TRUE, cc = TRUE, strata = NA){
  if (verbose){
    msg <- "Calculating ia for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  if (cc){
    rescc <- gid %>% clonecorrect(strata = strata) %>% poppr::ia()
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
  if (cc){
    ccout <- list(Iacc    = ~rescc["Ia"],
                  rbarDcc = ~rescc["rbarD"]
                  )
    out <- c(out, ccout)
  }
  if (keepdata){
    out <- c(out, list(dataset = ~list(gid)))
  }
  return(data_frame_(out))
}
