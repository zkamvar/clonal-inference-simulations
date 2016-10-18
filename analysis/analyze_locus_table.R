#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
"
================================================================================
Parse RDA files containing genclone objects, analyze locus diversity statistics, 
and save the resulting data frames.

Usage: analyze_and_save_ia.R [-dvm -o PATH] [FILE...]

Options:
 -h,--help              show this message and exit
 -v,--verbose           record progress
 -d,--debug             record EVERYTHING
 -m,--mutant            flag for indicating single mutant locus. If this is flagged, two analyses are run.
 -o PATH,--output=PATH  default path to place Rdata files [default: ~/rda_files]

================================================================================
" -> doc


# Processing command line arguments ---------------------------------------

opt <- docopt(doc, help = TRUE)

if (length(opt$FILE) == 0){
  message("\nFILE required. See options below\n")
  docopt(doc, help = TRUE, args = "-h")
}

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


# wrapper to exclude first locus before analyzing diversity.
tidy_no_first_locus <- function(x, ...){
  tidy_locus_table(x[loc = -1, mlg.reset = TRUE], ...)
}

# Processing data ---------------------------------------------------------

# For each file we will
# 1. Load the file into R (since it will be an Rda file)
# 2. Select the dataset column
# 3. apply the tidy_diversity function to each data set
# 4. store each result as a row in the table with pop as the key.
# 5. save the results to the file name, replacing data with "divtable".
for (f in opt$FILE){
  indat <- load(f) 
  if (opt$verbose) message(paste("Analyzing", nrow(get(indat)), "populations ... "))
  set.seed(opt$seed)
  res <- indat %>% 
    get() %>%
    select(dataset) %>% 
    apply(1, listfun, tidy_locus_table,  verbose = opt$debug) %>% 
    bind_rows()
  if (opt$mutant){
    if (opt$verbose) message(paste("Analyzing populations without first locus ... "))
    set.seed(opt$seed)
    mres <- indat %>% 
      get() %>%
      select(dataset) %>% 
      apply(1, listfun, tidy_no_first_locus,  verbose = opt$debug) %>% 
      bind_rows()
  }
  outf <- gsub(".DATA", ".locustable", basename(f))
  assign(x = outf, res)
  outf_location <- file.path(opt$output, outf)
  if (opt$verbose) message(paste("saving results to", outf_location))
  save(list = outf, file = outf_location)
  if (opt$mutant){
    moutf <- gsub(".DATA", ".locustable.mutant", basename(f))
    assign(x = moutf, mres)
    moutf_location <- file.path(opt$output, moutf)
    if (opt$verbose) message(paste("saving results to", moutf_location))
    save(list = moutf, file = moutf_location)
  }
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}

