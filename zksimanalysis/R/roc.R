#' Receiver Operator Characteristic
#'
#' @param df a data frame with p-values to compare
#' @param compare a two-element character vector specifying the levels of sexual reproduction to be the true positive and false positive.
#' @param stat What statistic should be compared
#' @param count.na When \code{TRUE}, missing values of \code{stat} are counted as positive
#' @param group the columns of data to group by
#' @param alpha the level at which to compare
#' @param flip a column of logicals used to flip a negative to a positive. This
#'   is used for extra criteria where the stat may not be considered significant
#'   based on alpha, but can be considered significant based on the presence of
#'   another factor
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' # single value
#' x <- roc(df = vals, group = c('sexrate'))
#'
#' # more granular
#' y <- roc(df = vals, group = c('sexrate', 'sample', 'mutation_rate'))
#'
#' # very granular
#' z <- roc(df = vals, group = c('sexrate', 'sample', 'mutation_rate', 'run', 'seed'))
#'
#' # Comparing clonal to sexual
#' res <- map_df(seq(0, 1, by = 0.01), roc, df = vals,
#'               group = c('sexrate', 'sample', 'mutation_rate'))
#' ggplot(res, aes(x = `False Positive`, y = `True Positive`, color = sample)) +
#'   geom_line() +
#'   geom_abline(slope = 1, lty = 2) +
#'   facet_wrap(~mutation_rate)
#' }
roc <- function(alpha = 0.05, df, compare = c("0.0000", "1.0000"),
                stat = "p.rD", augment = NULL, count.na = TRUE,
                group = c("sexrate", "run", "seed", "sample"),
                flip = NULL){
  flip <- if (is.null(flip)) "flip" else flip
  pf <- lazyeval::interp(~n_hits(stat, alpha, flip, count.na)/n(),
                         stat = as.name(stat),
                         flip = as.name(flip),
                         alpha = alpha,
                         count.na = count.na)
  score_column <- list(score = ~ifelse(sexrate == compare[1], "True Positive", "False Positive"))
  psd <- list(TPsd = ~roc_sd(`True Positive`, ntp),
              FPsd = ~roc_sd(`False Positive`, nfp))
  np  <- list(ntp  = ~sum(sexrate == compare[1]),
              nfp  = ~sum(sexrate == compare[2]),
              flip = FALSE)

  df %>%
    filter_(.dots = list(~sexrate %in% compare)) %>%
    mutate_(.dots = np) %>%
    group_by_(.dots = c(group, "ntp", "nfp")) %>%
    summarize_(.dots = list(positive_fraction = pf, alpha = alpha)) %>%
    mutate_(.dots = score_column) %>%
    ungroup() %>%
    mutate_(.dots = list(sexrate = ~compare[1])) %>%
    spread_(key_col = "score", value_col = "positive_fraction") %>%
    mutate_(.dots = psd) %>%
    select_(quote(-ntp), quote(-nfp))
}

n_hits <- function(stat, alpha, flip, count.na = TRUE){
  nas <- sum(is.na(stat)) # the flip will convert NaN to logical
  res <- (stat <= alpha)
  res <- sum(ifelse(!is.na(res) & alpha > 0, res | flip, res), na.rm = TRUE)
  res <- if (count.na) res + nas else res
  return(res)
}

roc_sd <- function(pf, n){
  sqrt((pf * (1 - pf))/(n - 1))
}
