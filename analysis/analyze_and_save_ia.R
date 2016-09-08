#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option(c("-f", "--files"), action = "store", type = "character",
                help = "Feather formatted file(s) to analyze"),
    make_option(c("-s", "--sample"), action = "store", type = "integer",
                default = c(10L, 25L, 50L, 100L),
                help = "Number of individuals to randomly sample from each population"),
    make_option(c("-p", "--permutations"), action = "store", type = "integer",
                default = 99L,
                help = "Number of permuatations to perform to calculate significance."),

)
