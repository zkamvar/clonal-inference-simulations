#' Generate random sample names
#'
#' Generates random sample names by sampling 10 letters with replacement.
#'
#' @param n the number of random names to generate
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
getNames <- function(n){
  vapply(1:n, function(x) paste(sample(letters, 10, replace = TRUE), collapse = ""), character(1))
}

#' Mutate a single SNP.
#'
#' One random SNP in one random chromosome. Used by sample_mutator.
#'
#' @param chrom a vector of raw elements
#' @param mutation a single integer representing a bit mask for a random
#'   mutation.
#'
#' @return a raw vector with one bit changed.
#'
#' @details the XOR gate is utilized to flip a single bit.
#'
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
snp_mutator <- function(chrom, mutation = sample(2^(0:7), 1)){
  posi        <- sample(length(chrom) - 1, 1) # sample chunk of 8 loci
  orig        <- as.integer(chrom[posi])      # convert to integer
  chrom[posi] <- as.raw(bitwXor(orig, mutation)) # Exclusive OR will flip the bit.
  return(chrom)
}

#' Mutate a single SNPbin object.
#'
#' propogate mutations randomly across chromosomes of a SNPbin object.
#'
#' @param snpbin a SNPbin object
#' @param mu the rate of mutation
#' @param nLoc the number of loci in the object.
#' @param rawchars an integer vector indicating the masking bits for mutation.
#'
#' @note The number of mutations is chosen from a single draw of a poisson
#'   distribution. Each mutation is propogated on a random chromosome in a
#'   random location.
#'
#' @return a mutated SNPbin object
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
sample_mutator <- function(snpbin, mu, nLoc, rawchars = 2^(0:7)){
  nmutations <- rpois(1, lambda = round(nLoc*mu))
  for (i in seq(nmutations)){
    chrom_index               <- sample(length(snpbin@snp), 1)
    chrom                     <- snpbin@snp[[chrom_index]]
    snpbin@snp[[chrom_index]] <- snp_mutator(chrom, sample(rawchars, 1))
  }
  return(snpbin)
}

#' Insert mutation into a genlight object
#'
#' This will randomly insert mutations at a given rate per sample.
#'
#' @param glt a genlight object
#' @param mu a mutation rate
#' @param samples a vector indicating which samples should be mutated
#'
#' @export
#' @return a mutated genlight object
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
pop_mutator <- function(glt, mu = 0.05, samples = TRUE){
  rawchrs <- 2^(0:7)
  glt@gen[samples] <- lapply(glt@gen[samples], sample_mutator, mu, nLoc(glt), rawchrs)
  return(glt)
}

#' Generate missing data in a single chromosome
#'
#' internal function utilized by NA_generator
#'
#' @param chrom a vector of type raw.
#' @param NA.posi the position of missing data across the chromosome. Note that
#'   these positions refer to positions of individual bits, not the individual
#'   raw elements.
#' @param rawchars a vector representing the masking bits. Passed for
#'   convenience.
#'
#' @details This is the trickiest part to generate missing data. Missing data in
#'   genlight objects are represented by an integer vector representing the
#'   missing position, but in the binary coding, all missing data are
#'   represented as "off" bits. When generating missing data, all bits in the
#'   missing positions must be flipped to the "off" position using an XOR gate.
#'   Given a single NA position, when challenged with an "on" bit in an OR gate,
#'   if the result is "on" the same as the input, then the bit must be flipped
#'   because the input was in the "on" position.
#'
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
NA_zeromancer <- function(chrom, NA.posi, rawchars = 2^(0:7)){
  nas <- ceiling(NA.posi/8) # Getting the location in the RAW vector
  zero_bits <- NA.posi %% 8 # Getting the location of the locus in a RAW element.
  zero_bits[zero_bits == 0] <- 8

  for (i in seq(length(nas))){
    eight_bit <- as.integer(chrom[nas[i]])
    the_locus <- rawchars[zero_bits[i]]
    # If the locus does not change in the OR, then there is a 1 and it needs to
    # be changed to a zero with XOR.
    if (eight_bit == bitwOr(eight_bit, the_locus)){
      chrom[nas[i]] <- as.raw(bitwXor(eight_bit, the_locus))
    }
  }
  return(chrom)
}

#' Generate missing data for a SNPbin object
#'
#' internal function utilized by pop_NA
#'
#' @param snpbin a SNPbin object
#' @param nloc the number of loci in the SNPbin object
#' @param na.perc the percentage of missing data to propogate
#' @param rawchars a vector representing masking bits. Passed for convenience.
#'
#' @details the number of missing loci is generated here by drawing from a
#'   poisson distribution.
#'
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
NA_generator <- function(snpbin, nloc, na.perc = 0.01, rawchars = 2^(0:7)){
  nas <- rpois(1, lambda = round(nloc*na.perc))
  NA.posi <- sort(sample(nloc, nas))
  for (i in seq(length(snpbin@snp))){
    snpbin@snp[[i]] <- NA_zeromancer(snpbin@snp[[i]], NA.posi, rawchars)
  }
  snpbin@NA.posi <- NA.posi
  return(snpbin)
}

#' Propogate missing data in a genlight object
#'
#' @param na.perc the percentage of missing data per sample
#' @param parallel a logical specifying whether or not parallel processing
#'   should be utilized
#' @param n.cores the number of cores to use with parallel processing.
#'
#' @return a genlight object with missing data propogated.
#'
#' @export
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
pop_NA <- function(glt, na.perc = 0.01, parallel = require('parallel'), n.cores = 2L){
  rawchars <- 2^(0:7)
  if (parallel){
    glt@gen <- mclapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars,
                        mc.cores = getOption("mc.cores", n.cores))
  } else {
    glt@gen <- lapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars)
  }

  return(glt)
}

#' Create crossover in a SNPbin object
#'
#' @param snpbin a SNPbin object
#'
#' @return a SNPbin object with its chromosomes crossed over
#'
#' @details this will randomly choose a cut point in the chromosomes of a
#'   diploid snpbin object and create a crossover event.
#'
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
crossover <- function(snpbin){
  chr1 <- snpbin@snp[[1]]
  chr2 <- snpbin@snp[[2]]
  chrlen <- length(chr1)
  cutpoint <- sample(2:(chrlen - 1), 1)
  first <- 1:(cutpoint - 1)
  last  <- cutpoint:chrlen
  snpbin@snp[[1]] <- c(chr1[first], chr2[last])
  snpbin@snp[[2]] <- c(chr2[first], chr1[last])
  return(snpbin)
}

#' mate two individuals from a genlight object
#'
#' @param snpbin a genlight object (yes, I know its mis-named).
#' @param ind1 the index of parent 1
#' @param ind2 the index of parent 2
#' @return a SNPbin object
#'
#' @details this will take two parental samples (they do not have to be
#'   different), have them each undergo crossover, and then randomly sample one
#'   chromosome from each to create the offspring.
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
mate <- function(snpbin, ind1, ind2){
  snpbin@gen[[ind1]] <- crossover(snpbin@gen[[ind1]])
  snpbin@gen[[ind2]] <- crossover(snpbin@gen[[ind2]])
  snpout <- snpbin@gen[[ind1]]
  snpout@snp[[1]] <- snpbin@gen[[ind1]]@snp[[sample(1:2, 1)]]
  snpout@snp[[2]] <- snpbin@gen[[ind2]]@snp[[sample(1:2, 1)]]
  return(snpout)
}

#' randomly mate a genlight object once
#'
#' @param glt a genlight object
#' @param err the number of mutations per mating event.
#'
#' @details Samples are chosen with replacement from the population twice to
#'   create each set of parents. They are then mated and mutated. See getSims
#'   for details
#'
#' @return a genlight object with the same number of individuals as glt
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
#' @keywords internal
random_mate <- function(glt, err){
  mat_pair <- matrix(integer(1), nrow = nInd(glt), ncol = 2)
  mat_pair[, 1] <- sample(nInd(glt), replace = TRUE)
  mat_pair[, 2] <- sample(nInd(glt), replace = TRUE)
  res <- apply(mat_pair, 1, function(x) mate(glt, x[1], x[2]))
  glt@gen <- res
  if (err > 0) glt <- pop_mutator(glt, err)
  return(glt)
}

#' randomly mate a genlight object for n generations.
#'
#' @param glt a genlight object
#' @param err the number of mutations per mating event
#' @param gen the number of generations to mate the same population
#'
#' @return a genlight object
#' @author Zhian N. Kamvar
#' @references ZN Kamvar, JC Brooks, and NJ Grünwald. 2015. Supplementary
#'   Material for Frontiers Plant Genetics and Genomics 'Novel R tools for
#'   analysis of genome-wide population genetic data with emphasis on
#'   clonality'. DOI:
#'   \href{http://dx.doi.org/10.5281/zenodo.17424}{10.5281/zenodo.17424}
#'
#'   Kamvar ZN, Brooks JC and Grünwald NJ (2015) Novel R tools for analysis of
#'   genome-wide population genetic data with emphasis on clonality. Front.
#'   Genet. 6:208. doi:
#'   \href{http://dx.doi.org/10.3389/fgene.2015.00208}{10.3389/fgene.2015.00208}
random_mate_gen <- function(glt, err = 5e-3, gen = 1){
  for (i in seq(gen)){
    glt <- random_mate(glt, err)
  }
  return(glt)
}
#' Helper function to display bitcode.
#'
#' @param y a raw string
#' @return a character string of 1s and 0s
#' @noRd
#' @example
#' # Bitshifters
#' shifts <- as.raw(2^(0:7))
#' binary_char_from_hex(shifts)
binary_char_from_hex <- function(y){
  vapply(y, function(x) paste(as.integer(rawToBits(x)), collapse=""), character(1))
}
