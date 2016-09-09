#!/usr/bin/env Rscript
#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("zksimanalysis"))

option_list <- list(
    make_option(c("-f", "--files"), action = "store", type = "character",
                default = NULL,
                help = "Feather formatted file(s) to analyze"),
    make_option(c("-s", "--seed"), action = "store", type = "integer",
                default = 20160909),
    make_option(c("-n", "--nsample"), action = "store", type = "integer",
                default = c(10L, 25L, 50L, 100L),
                help = "Number of individuals to randomly sample from each population"),
    make_option(c("-p", "--permutations"), action = "store", type = "integer",
                default = 99L,
                help = "Number of permuatations to perform to calculate significance."),
    make_option(c("-l", "--locus"), action = "store", type = "character",
                default = "Locus",
                help = "Prefix for identifying locus columns")
)

opt_parser <- OptionParser(option_list = option_list)
opt        <- parse_args(opt_parser)

if (is.null(opt$files)){
  stop("must specify a file")
}

for (f in opt$files){
  set.seed(opt$seed)
  datalist <- lapply(f, feather2genind, locus_prefix = opt$locus, sample = opt$nsample) %>%
    lapply(setPop, ~pop) %>%
    lapply(seppop) %>%
    lapply(lapply(tidy_ia, sample = opt$permutations, hist = FALSE, quiet = TRUE) %>% bind_rows()) %>%
    bind_rows()
}
