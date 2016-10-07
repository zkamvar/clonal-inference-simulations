#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
"
================================================================================
Parse RDA files containing genclone objects, analyze locus contribution to
observed multilocus genotypes. This uses the analyze_locus_contribution()
function from zksimanalysis and will calculate: A genotype accumulation curve,
round-robin multilocus genotypes, diversity statistics based on the number of
multilocus genotypes contributed by each locus, and probability of genotypes for
each sample.

Usage: analyze_locus_contribution.R [-dv -s SEED -p PERMUTATIONS -o PATH] [FILE...]

Options:
 -h,--help                                    show this message and exit
 -v,--verbose                                 record progress
 -d,--debug                                   record EVERYTHING
 -s SEED,--seed=SEED                          random seed [default: 20160909]
 -p PERMUTATIONS,--permutations=PERMUTATIONS  number of permutations for the index of association [default: 99]
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

if (!dir.exists(opt$output)){
  dir.create(opt$output)



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
# Number of multilocus genotypes for a matrix of genotypes
# Number of multilocus genotypes for a matrix of genotypes
NMLG <- function(x){
  if (length(dim(x)) > 1){
    N <- rep(ncol(x), nrow(x))
  } else {
    N <- length(x)
  }
  return(N)
}

# Unbiased simpson's diversity is corrected for sample size.
uSimp <- function(x){
  lambda <- vegan::diversity(x, "simpson")
  N      <- NMLG(x)
  return((N/(N-1))*lambda)
}

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
    apply(1, listfun, tidy_locus_contribution, 
          uSimp = uSimp, 
          n = opt$permutations, 
          verbose = opt$debug) %>% 
    bind_rows()
  outf <- gsub(".DATA", ".contrib", basename(f))
  assign(x = outf, res)
  outf_location <- file.path(opt$output, outf)
  if (opt$verbose) message(paste("saving results to", outf_location))
  save(list = outf, file = outf_location)
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}

