#' Generate a tidy data frame of the index of association
#'
#' @param gid a genind or genclone object
#' @param ... any parameters to be passed to \code{\link[poppr]{ia}}
#' @param verbose when \code{TRUE}, a status message will be printed at the
#'   beginning of the analysis. \code{FALSE} suppresses this message.
#' @param keepdata when \code{TRUE}, the data set is kept in a column called
#'   "dataset"
#' @param cc return values of the clone-corrected data set.
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
tidy_ia <- function(gid, ..., verbose = TRUE, keepdata = TRUE, cc = TRUE,
                    strata = NA){
  if (verbose){
    msg <- "Calculating ia for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  if (inherits(gid, "genlight")){
    return(tidy_genomic_ia(gid, ..., keepdata = keepdata, cc = cc, strata = strata))
  }
  if (cc){
    rescc <- gid %>%
      clonecorrect(strata = strata) %>%
      poppr::ia(..., valuereturn = TRUE)
  }
  res  <- poppr::ia(gid, ..., valuereturn = TRUE)
  pops <- popNames(gid)
  pops <- if (length(pops) > 1) list(pops) else pops
  out  <- process_ialist(res)
  if (cc){
    ccout        <- process_ialist(rescc)
    names(ccout) <- paste0(names(ccout), "cc")
    out          <- c(out, ccout)
  }
  if (keepdata){
    out <- c(out, list(dataset = ~list(gid)))
  }
  out <- c(out, list(pop = ~pops))
  return(data_frame_(out))
}

#' Generate tidy data frame of jackknife statistics of the index of association
#'
#' @inheritParams tidy_ia
#' @param ... parameters to be passed on to \code{\link[poppr]{jack.ia}}
#'
#' @return a tibble
#' @export
#'
#' @examples
#' library('poppr')
#' data(partial_clone)
#' strata(partial_clone) <- data.frame(pop = pop(partial_clone))
#' partial_clone %>%
#'   seppop() %>%
#'   lapply(tidy_jackia, reps = 999, strata = 1) %>%
#'   bind_rows()
tidy_jackia <- function(gid, ..., verbose = TRUE, strata = NA){
  if (verbose){
    msg <- "Calculating jackknife stats for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  res <- vector(mode = "list", length = 2) %>% setNames(c("observed", "samples"))
  obs <- gid %>%
    poppr::clonecorrect(strata = strata)%>%
    poppr::ia()
  res[["samples"]] <- try(poppr::jack.ia(gid, ...))
  if (inherits(res[["samples"]], "try-error")){
    out <- list(
      jack.p.Ia = ~NA_real_,
      jack.p.rD = ~NA_real_,
      jack.ia   = ~list(NA_real_),
      jack.rd   = ~list(NA_real_),
      pop       = ~NA_character_
    )
    return(data_frame_(out))
  }
  res[["observed"]] <- c(obs["Ia"],
                         p.Ia = makep(obs["Ia"], res[["samples"]][["Ia"]]),
                         obs["rbarD"],
                         p.rD = makep(obs["rbarD"], res[["samples"]][["rbarD"]])
  )
  class(res) <- "ialist"
  pops <- popNames(gid)
  pops <- if (length(pops) > 1) list(pops) else pops
  out <- process_ialist(res)[-c(1, 3)]
  out <- c(out, list(pop = ~pops))
  data_frame_(out) %>%
    setNames(c("jack.p.Ia", "jack.p.rD", "jack.ia", "jack.rd", "pop"))
}

# helper function to perform a one-sided test
makep <- function(x, y) (sum(x >= y, na.rm = TRUE) + 1)/(length(y) + 1)

process_ialist <- function(res){
  if (inherits(res, "ialist")){
    vals <- res[[1]]
    out  <- list(Ia         = ~vals["Ia"],
                 p.Ia       = ~vals["p.Ia"],
                 rbarD      = ~vals["rbarD"],
                 p.rD       = ~vals["p.rD"],
                 samples.ia = ~list(res[[2]]$Ia),
                 samples.rd = ~list(res[[2]]$rbarD)
                 )
  } else {
    out <- list(Ia    = ~res["Ia"],
                rbarD = ~res["rbarD"]
                )
  }
  return(out)
}

tidy_genomic_ia <- function(gid, ..., keepdata, cc, strata){
  res  <- genomic_ia(gid, ...)
  pops <- popNames(gid)
  pops <- if (length(pops) > 1) list(pops) else pops
  if (inherits(res, "list")){
    out  <- list(rbarD   = ~res$observed["rbarD"],
                 p.rD    = ~res$observed["p.rD"],
                 samples = ~list(res$samples)
                 )
  } else {
    out  <- list(rbarD = ~res["rbarD"])
  }

  if (keepdata){
    out <- c(out, list(dataset = ~list(gid)))
  }
  out <- c(out, list(pop = ~pops))
  return(data_frame_(out))
}
