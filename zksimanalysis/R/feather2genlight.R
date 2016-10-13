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
feather2genlight <- function(ff, locus_regex = "^[0-9]+?_[0-9]+?$", sample = 50, snpclone = TRUE, verbose = FALSE){
  if (verbose) message(paste0("Reading ", ff, " ..."))
  if (!is.null(sample)){
    fdf  <- feather::read_feather(ff) %>% group_by_(~pop) %>% sample_n(sample)
  } else {
    fdf <- feather::read_feather(ff) %>% group_by_(~pop)
  }
  if (verbose) message(paste0("Converting to genlight/snpclone ..."))
  loci <- fdf %>% ungroup() %>% select_(~matches(locus_regex))
  sta  <- fdf %>% select_(~-matches(locus_regex))
  the_loci  <- names(loci)
  chrom_loc <- vapply(the_loci, function(i) strsplit(i, "_")[[1]], character(2))
  res <- new("genlight", loci,
             chromosome = chrom_loc[1, ],
             position = get_position(chrom_loc[2, ]),
             strata = sta
             )
  res <- if (snpclone) as.snpclone(res)
  if (verbose) message("Done.")
  return(res)
}


get_position <- function(x){
  charlen  <- nchar(x)
  ismax    <- charlen == max(charlen)
  x[ismax] <- substr(x[ismax], 2, max(charlen))
  return(as.integer(x))
}
