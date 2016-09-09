#' Convert a feather file to a genind object
#'
#' @param ff a feather file on disk
#' @param locus_prefix An identifier for the locus columns (Default: "Locus")
#' @param genclone convert to a genclone object (Default: TRUE)
#' @param verbose report progress of import? (Default: FALSE)
#' @param sample the number of individuals to sample per population. (Default: 50)
#'
#' @return a genind object
#' @export
#'
#' @examples
#' library("tidyr")
#' ex  <- example_data()
#' gid <- feather2genind(ex)
#' # split the population name into its components
#' ex_run  <- "twenty_loci([0-9]+?)/"
#' ex_seed <- "seed_([0-9]+?)_"
#' ex_sex  <- "sex_([0-9.]+?)_"
#' ex_gen  <- "gen_([0-9]+?)_"
#' ex_rep  <- "rep_([0-9]+?).pop"
#' strata(gid) <- strata(gid) %>%
#'   extract(pop, c("run", "seed", "sexrate", "gen", "rep"),
#'           paste0(ex_run, ex_seed, ex_sex, ex_gen, ex_rep),
#'           remove = FALSE)
#' setPop(gid) <- ~run/rep/seed
#' gid
feather2genind <- function(ff, locus_prefix = "Locus", sample = 50, genclone = TRUE, verbose = FALSE){
  if (verbose) message(paste0("Reading ", ff, " ..."))
  fdf  <- feather::read_feather(ff) %>% group_by_(~pop) %>% sample_n(sample)
  if (verbose) message(paste0("Converting to genind/genclone ..."))
  loci <- fdf %>% select_(~starts_with(locus_prefix))
  sta  <- fdf %>% select_(~-starts_with(locus_prefix))
  gid  <- df2genind(loci, sep = "/", strata = sta)
  gid  <- if (genclone) as.genclone(gid)
  if (verbose) message("Done.")
  return(gid)
}


#' permute the index of assocation and return a tibble
#'
#' @param gid a genind or genclone object
#' @param ... any parameters to be passed to \code{\link[poppr]{ia}}
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
tidy_ia <- function(gid, ..., verbose = TRUE){
  if (verbose){
    msg <- "Calculating ia for"
    if (nPop(gid) == 1){
      msg <- paste(msg, popNames(gid), "...")
    } else {
      msg <- paste(msg, nInd(gid), "samples ...")
    }
    message(msg)
  }
  res <- poppr::ia(gid, ..., valuereturn = TRUE)
  vals <- res[[1]]
  pops <- popNames(gid)
  pops <- if (length(pops) > 1) list(pops) else pops
  out <- list(Ia         = ~vals["Ia"],
              p.Ia       = ~vals["p.Ia"],
              rbarD      = ~vals["rbarD"],
              p.rD       = ~vals["p.rD"],
              samples.ia = ~list(res[[2]]$Ia),
              samples.rd = ~list(res[[2]]$rbarD),
              pop        = ~pops
              )
  return(data_frame_(out))
}

#' A wrapper to get example data
#'
#' @return a string with the path to the example
#' @export
#'
#' @examples
#' example_data()
example_data <- function(){
  system.file("files", "sample_data.feather", package = "zksimanalysis")
}
