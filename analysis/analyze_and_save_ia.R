#!/usr/bin/env Rscript
#
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("zksimanalysis"))
suppressPackageStartupMessages(library("purrr"))
"Usage: analyze_and_save_ia.R [-s SEED [-n NSAMPLE...] -p PERMUTATIONS -l LOCUS] -o PATH -f FILES...

Options:
 -s SEED, --seed=SEED
    random seed [default: 20160909]
 -n NSAMPLE..., --nsample=NSAMPLE...
    number of samples/population [default: 10 25 50 100]
 -p PERMUTATIONS, --permutations=PERMUTATIONS
    number of permutations for the index of association [default: 99]
 -l LOCUS, --locus=LOCUS
    prefix to identify the locus columns [default: Locus]
 -o PATH, --outdir=PATH
    the directory in which to store the rdata files [default: ~/rda_files]
 -f FILES..., --files=FILES...
    feather-format files to read" -> doc

docopt(doc)

if (is.null(opt$files)){
  stop("must specify a file")
}

totsamp <- sum(opt$nsample)
if (!dir.exists(opt$outdir)){
  dir.create(opt$outdir)
}

for (f in opt$files){
  set.seed(opt$seed)
  indat <- feather2genind(f, locus_prefix = opt$locus, sample = totsamp) %>%
    setPop(~pop)
  subpops <- rep(opt$nsample, nPop(indat)) %>%
    rep(., .) %>%
    paste0("sam_", .)
  res <- indat %>%
    addStrata(data.frame(sample_size = subpops)) %>%
    setPop(~pop/sample_size) %>%
    seppop() %>%
    map(tidy_ia, sample = opt$permutations, hist = FALSE, quiet = FALSE) %>%
    bind_rows()
  outf <- make.names(f)
  assign(x = outf, res)
  outf_location <- paste0(opt$outdir, "/", outf, ".rda")
  save(list = outf, file = outf_location)
}
