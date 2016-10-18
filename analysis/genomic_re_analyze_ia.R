#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
suppressPackageStartupMessages(library("stringr"))
"
================================================================================
Parse feather formatted data, analyze the index of association on the whole data
set, propogate missing data, analyze the index of association there, and save
the result.

Usage: genomic_analyze_and_save_ia.R [-hdvk [-m MISSING...] -b BLOCKSIZE -s SEED [-n NSAMPLE...] -p PERMUTATIONS -l LOCUS -o PATH] [FILE...]

Options:
 -h,--help                                    show this message and exit
 -v,--verbose                                 record progress
 -d,--debug                                   record EVERYTHING
 -m MISSING...,--missing=MISSING...           include analyze missing data percentage (no permutation) [default: 0.01 0.05 0.10]
 -s SEED,--seed=SEED                          random seed [default: 20160909]
 -b BLOCK,--blocksize=BLOCK                   number of blocks by which to sample [default: 1]
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
opt$blocksize    <- as.integer(opt$blocksize)
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
suppressPackageStartupMessages(library("tidyverse"))

# Processing data ---------------------------------------------------------

for (f in opt$FILE){
  set.seed(opt$seed) # same seed for each file
  indat <- load(opt$FILE)
  if (opt$verbose) message(paste("Analyzing", indat, "... "))
   res <- indat %>%
    get() %>%
    select(dataset) %>%
    apply(1, listfun, tidy_ia,
          sample    = opt$permutations,
          keepdata  = FALSE,
          blocksize = opt$blocksize,
          verbose   = opt$verbose,
          quiet     = !opt$debug) %>%
    bind_rows()
  # Saving everything.
  outf <- gsub(".DATA", ".shufflechrom", basename(f))
  assign(x = outf, res)
  outf_location <- file.path(opt$output, outf)
  if (opt$verbose) message(paste("saving results to", outf_location))
  save(list = outf, file = outf_location)
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}
