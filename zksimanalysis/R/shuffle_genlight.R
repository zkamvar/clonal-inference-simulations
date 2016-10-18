#' resample a genlight/genclone object and produce a
#'
#' @param x a snpclone or genlight object
#' @param mat a matrix derived from x
#' @param blocksize the size of the blocks at which genotypes should be shuffled.
#'
#' @return a snpclone or genlight object
#' @export
#'
#' @examples
#' gl <- glSim(100, 100, 100, parallel = FALSE)
#' gl
#' bitwise.ia(gl)
#' gls <- shuffle_genlight(gl, as.matrix(gl))
#' gls
#' bitwise.ia(gls)
shuffle_genlight <- function(mat, x, blocksize = 100){

  if (blocksize == 1){
    shuff <- new("genlight", apply(mat, 2, sample), parallel = FALSE)
  } else {
    if (is.null(position(x))){
      position(x) <- seq(nLoc(x))
    }
    pos <- position(x)
    if (!is.null(chromosome(x))){
      # Each chromosome has it's own position.
      # In this case, get large round number for each chromosome break.
      maxp <- 10^ceiling(log(max(pos), 10))
      if (length(unique(pmin(pos))) < nLoc(x)){
        pos <- split(pos, chromosome(x))
        for (p in seq(pos)){
          if (p != 1){
            pos[[p]] <- pos[[p]] + (maxp*(p - 1))
          }
        }
        pos <- unlist(pos, use.names = FALSE)
      }
    }
    blocks   <- seq(blocksize, max(pos), by = blocksize)
    blockmat <- matrix(blocks, nrow = length(blocks), ncol = 2)
    blockmat[, 1] <- blockmat[, 1] - blocksize + 1

    for (i in seq(blocks)){
      blockrange <- pos %in% blockmat[i, 1]:blockmat[i, 2]
      if (any(blockrange)){
        s <- sample(nrow(mat))
        mat[, blockrange] <- mat[s, blockrange]
      }
    }
  }
  shuff <- new("genlight", mat, parallel = FALSE)
  x@gen <- shuff@gen
  return(x)
}

#' Calculate index of association from reshuffled genlight object
#'
#' @param x a genlight or snpclone object
#' @param mat a matrix derived from x
#' @param blocksize the size of the blocks at which genotypes should be shuffled.
#'
#' @return a value of the index of association
#' @export
#'
#' @examples
#' gl <- glSim(100, 100, 100, parallel = FALSE)
#' bitwise.ia(gl)
#' sample_bitwise_ia(gl, as.matrix(gl))
sample_bitwise_ia <- function(mat, x, blocksize = 100){
  bitwise.ia(shuffle_genlight(mat, x, blocksize), threads = 1L)
}

#' A test for significance with genomic ia
#'
#' @param x a genlight or snpclone object
#' @param sample the number of samples for significance testing
#' @param blocksize the size of the blocks at which genotypes should be shuffled.
#'
#' @return a list with the observed value and pvalue and the samples used.
#' @export
#'
#' @note This is a slow procedure. When tested on 185 samples over 10,000 SNPs,
#' this took about 5 seconds/sample, which is 84 minutes of computing time on
#' my mac.
#'
#' @examples
#' gl <- glSim(100, 200, 0, parallel = FALSE)
#' genomic_ia(gl, sample = 999)
genomic_ia <- function(x, sample = 0, blocksize = 100, quiet = FALSE){
  obs  <- bitwise.ia(x, threads = 1L)
  mat  <- as.matrix(x)
  if (sample == 0) return(c(rbarD = obs))
  exp  <- vector(mode = "numeric", length = sample)
  if (!quiet) p <- dplyr::progress_estimated(sample)
  for (i in seq(sample)){
    exp[i] <- sample_bitwise_ia(mat, x, blocksize)
    if (!quiet) p$tick()$print()
  }
  if (!quiet) cat("\n")
  pval <- (sum(exp >= obs) + 1)/(sample + 1)
  return(list(observed = c(rbarD = obs, p.rD = pval), samples = exp))
}
