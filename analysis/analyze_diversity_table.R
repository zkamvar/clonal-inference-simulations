#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("docopt"))
"
================================================================================
Parse RDA files containing genclone objects, analyze diversity statistics, and
save the resulting data frames.

Usage: analyze_and_save_ia.R [-dv -s SEED -p PERMUTATIONS -o PATH] [FILE...]

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
# Number of multilocus genotypes for a matrix of genotypes
NMLG <- function(x){
  if (length(dim(x)) > 1){
    N <- rowSums(x > 0)
  } else {
    N <- sum(x > 0)
  }
  return(N)
}

# Unbiased simpson's diversity is corrected for sample size.
uSimp <- function(x){
  lambda <- vegan::diversity(x, "simpson")
  N      <- NMLG(x)
  return((N/(N-1))*lambda)
}

# This estimates the slope of a log-log power law linear model.
power_law_beta <- function(x){
  if (length(x) == 1) return(NA_real_)
  xpow <- displ(x[x > 0])                  # Generate the distribution
  xpow$setPars(estimate_pars(xpow))        # Estimate the parameters
  xdat <- plot(xpow, draw = FALSE)         # Extract the data
  xlm  <- lm(log(y) ~ log(x), data = xdat) # Run log-log linear model for slope
  return(-coef(xlm)[2])
}

# Wrapper for the estimate.
Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
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
    apply(1, listfun, tidy_diversity, 
          B = Beta, 
          uSimp = uSimp, 
          NMLG = NMLG,
          n = opt$permutations, 
          verbose = opt$debug) %>% 
    bind_rows()
  outf <- gsub(".DATA", ".divtable", basename(f))
  assign(x = outf, res)
  outf_location <- file.path(opt$output, outf)
  if (opt$verbose) message(paste("saving results to", outf_location))
  save(list = outf, file = outf_location)
}
if (opt$verbose){
  options(width = 100)
  print(devtools::session_info())
}

