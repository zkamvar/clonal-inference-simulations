#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
suppressPackageStartupMessages(library("stringr"))
"
================================================================================
Parse feather formatted data, analyze the index of association and save the result

Usage: analyze_and_save_ia.R [-hdvk [-m MISSING...] -s SEED [-n NSAMPLE...] -p PERMUTATIONS -l LOCUS -o PATH] [FILE...]

Options:
 -h,--help                                    show this message and exit
 -v,--verbose                                 record progress
 -d,--debug                                   record EVERYTHING
 -k,--keep                                    keep the data in the resulting data frame
 -m MISSING...,--missing=MISSING...           include analyze missing data percentage (no permutation) [default: 0.01 0.05 0.10]
 -s SEED,--seed=SEED                          random seed [default: 20160909]
 -n NSAMPLE...,--nsample=NSAMPLE...           number of samples/population (single quoted list of integers) [default: 10 25 50 100]
 -p PERMUTATIONS,--permutations=PERMUTATIONS  number of permutations for the index of association [default: 99]
 -l LOCUS,--locus=LOCUS                       prefix to identify the locus columns [default: ^[0-9]+?_[0-9]+?$]
 -o PATH,--output=PATH                        default path to place Rdata files [default: ~/genomic_rda_files]

================================================================================
" -> doc


# Processing command line arguments ---------------------------------------

opt <- docopt(doc, help = TRUE)
print(opt)

if (length(opt$FILE) == 0){
  message("\nFILE required. See options below\n")
  docopt(doc, help = TRUE, args = "-h")
}

opt$seed         <- as.integer(opt$seed)
opt$permutations <- as.integer(opt$permutations)
opt$nsample      <- as.integer(str_split(opt$nsample, "[[:blank:]]")[[1]])
print(opt$missing)
opt$missing      <- as.numeric(str_split(opt$missing, "[[:blank:]]")[[1]])

totsamp <- sum(opt$nsample)
if (!dir.exists(opt$output)){
  dir.create(opt$output)
}


# Initializing packages ---------------------------------------------------

if (opt$verbose){
  the_options <- capture.output(print(opt))
  the_options <- the_options[grepl("^ \\$ [[:alnum:]].+?$", the_options)]
  message("The options:\n\n")
  message(paste0(the_options, "\n"))
  message("\nLoading packages ...\n")
}
suppressPackageStartupMessages(library("zksimanalysis"))
suppressPackageStartupMessages(library("purrr"))



# Helper Functions --------------------------------------------------------

# Full analysis here takes the data, separates, populations, and analyzes
# $\bar{r}_d$ with permutations
full_analysis <- . %>%
  seppop() %>%
  map(tidy_ia,
      sample   = opt$permutations,
      keepdata = opt$keep,
      verbose  = opt$verbose,
      quiet    = !opt$debug) %>%
  bind_rows()

# Missing analysis will calculate the value of $\bar{r}_d$ for all the data
# sets with missing data and bind them into a tibble.
missing_analysis_match <- . %>%
  map(seppop, parallel = FALSE) %>%
  map(map_dbl,
      bitwise.ia,
      threads = 1L) %>%
  bind_rows()

missing_analysis_nomatch <- . %>%
  map(seppop, parallel = FALSE) %>%
  map(map_dbl,
      bitwise.ia,
      missing_match = FALSE,
      threads = 1L) %>%
  bind_rows()

# Processing data ---------------------------------------------------------

for (f in opt$FILE){
  set.seed(opt$seed) # same seed for each file
  indat <- feather2genlight(f, locus_regex = opt$locus, sample = totsamp, verbose = opt$verbose) %>%
    setPop(~pop)
  inpops  <- nPop(indat)
  allpops <- inpops * length(opt$nsample)
  subpops <- rep(opt$nsample, inpops) %>%
    rep(., .) %>%
    paste0("sam_", .)

  # Processing the incoming data so that it can be split into subpopulations.
  indat <- indat %>%
    addStrata(data.frame(sample_size = subpops)) %>%
    setPop(~pop/sample_size)

  # Set missing will create a list of data sets with missing data. The input is
  # the vector of missing percentages.
  set_missing <- . %>% {
    set.seed(opt$seed)
    pop_NA(indat, na.perc = ., parallel = FALSE)
  }
  miss <- map(opt$missing, set_missing) %>%
    setNames(paste0("missing_", opt$missing))

  if (opt$verbose) message(paste("Analyzing", allpops, "populations ... "))
  # Analyze $\bar{r}_d$ with permutations
  res  <- indat %>%
    full_analysis
  if (opt$verbose) message(paste("Analyzing missing data (with matching) ... "))
  # Analyze $\bar{r}_d$ with missing data (missing data always matches)
  mres <- miss %>%
    missing_analysis_match
  if (opt$verbose) message(paste("Analyzing missing data (no matching) ... "))
  # Analyze $\bar{r}_d$ with missing data (missing data never matches)
  mres_nomatch <- miss %>%
    setNames(paste0("missing_nomatch_", opt$missing)) %>%
    missing_analysis_nomatch
  # Preserve the augmented data sets.
  mdat <- miss %>%
    setNames(paste0("dataset_", opt$missing)) %>%
    map(seppop, parallel = FALSE) %>%
    bind_rows()
  res <- bind_cols(mres, mres_nomatch, res, mdat)

  # Saving everything.
  outf <- make.names(f)
  if (opt$keep){
    outfDATA <- paste0(outf, ".DATA")
    resDATA  <- res %>% dplyr::select(matches("dataset"), matches("pop"))
    assign(x = outfDATA, resDATA)
    outf_data_location <- paste0(opt$output, "/", outfDATA, ".rda")
    if (opt$verbose) message(paste("saving data to", outf_data_location))
    save(list = outfDATA, file = outf_data_location)
  }
  res  <- res %>% dplyr::select(-matches("dataset"))
  assign(x = outf, res)
  outf_location <- paste0(opt$output, "/", outf, ".rda")
  if (opt$verbose) message(paste("saving results to", outf_location))
  save(list = outf, file = outf_location)
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}
