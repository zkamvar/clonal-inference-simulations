#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
suppressPackageStartupMessages(library("stringr"))
"
================================================================================
Parse feather formatted data, analyze the index of association and save the result

Usage: analyze_and_save_ia.R [-dvkcm -s SEED [-n NSAMPLE...] -p PERMUTATIONS -l LOCUS -o PATH] [FILE...]

Options:
 -h,--help                                    show this message and exit
 -v,--verbose                                 record progress
 -d,--debug                                   record EVERYTHING
 -m,--mutant                                  flag for indicating single mutant locus. If this is flagged, two analyses are run.
 -k,--keep                                    keep the data in the resulting data frame
 -c,--clonecorrect                            add values from clone-correction
 -s SEED,--seed=SEED                          random seed [default: 20160909]
 -n NSAMPLE...,--nsample=NSAMPLE...           number of samples/population (single quoted list of integers) [default: 10 25 50 100]
 -p PERMUTATIONS,--permutations=PERMUTATIONS  number of permutations for the index of association [default: 99]
 -l LOCUS,--locus=LOCUS                       prefix to identify the locus columns [default: Locus]
 -o PATH,--output=PATH                        default path to place Rdata files [default: ~/rda_files]

================================================================================
" -> doc


# Processing command line arguments ---------------------------------------

opt <- docopt(doc, help = TRUE)

if (length(opt$FILE) == 0){
  message("\nFILE required. See options below\n")
  docopt(doc, help = TRUE, args = "-h")
}

opt$seed         <- as.integer(opt$seed)
opt$permutations <- as.integer(opt$permutations)
opt$nsample      <- as.integer(str_split(opt$nsample, "[[:blank:]]")[[1]])

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


# Processing data ---------------------------------------------------------
#
# Here, each file is processed by:
#  1. sampling a sum of all the samples from the file
#     This is done as a sum of all the sample sizes to reduce read time and
#     ensure sampling without replacement.
#  2. adding a strata for each sample size
#  3. splitting the populations by the sample sizes
#  4. calculating the index of association with permutation for each population
#  5. binding these results into a tidy data frame
#  6. saving this as a file in the output directory
for (f in opt$FILE){
  set.seed(opt$seed) # same seed for each file
  indat <- feather2genind(f, locus_prefix = opt$locus, sample = totsamp, verbose = opt$verbose) %>%
    setPop(~pop)
  inpops  <- nPop(indat)
  allpops <- inpops * length(opt$nsample)
  subpops <- rep(opt$nsample, inpops) %>%
    rep(., .) %>%
    paste0("sam_", .)
  if (opt$verbose) message(paste("Analyzing", allpops, "populations ... "))
  res <- indat %>%
    addStrata(data.frame(sample_size = subpops)) %>%
    setPop(~pop/sample_size) %>%
    seppop() %>%
    map(tidy_ia,
        sample   = opt$permutations,
        keepdata = opt$keep,
        cc       = opt$clonecorrect,
        strata   = ~pop/sample_size,
        verbose  = opt$verbose,
        hist     = FALSE,
        quiet    = !opt$debug) %>%
    bind_rows()
  if (opt$mutant){
  # At this point, if the first locus is a locus with a higher mutation rate
  # than the other loci, It is removed, MLGS recalculated, and everything is
  # processed again.
    if (opt$verbose) message(paste("Analyzing populations without first locus ... "))
    mres <- indat[loc = -1, mlg.reset = TRUE] %>%
      addStrata(data.frame(sample_size = subpops)) %>%
      setPop(~pop/sample_size) %>%
      seppop() %>%
      map(tidy_ia,
          sample   = opt$permutations,
          keepdata = FALSE,
          cc       = opt$clonecorrect,
          strata   = ~pop/sample_size,
          verbose  = opt$verbose,
          hist     = FALSE,
          quiet    = !opt$debug) %>%
      bind_rows()
  }
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

  if (opt$mutant){
    moutf <- paste0(outf, ".mutant")
    assign(x = moutf, mres)
    moutf_location <- paste0(opt$output, "/", moutf, ".rda")
    if (opt$verbose) message(paste("saving mutant results to", moutf_location))
    save(list = moutf, file = moutf_location)
  }
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}
