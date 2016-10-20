#' Receiver Operator Characteristic
#'
#' @param df a data frame with p-values to compare
#' @param compare a two-element character vector specifying the levels of sexual reproduction to be the true positive and false positive.
#' @param stat What statistic should be compared
#' @param count.na When \code{TRUE}, missing values of \code{stat} are counted as positive
#' @param group the columns of data to group by
#' @param alpha the level at which to compare
#'
#' @return a data frame
#' @export
#'
#' @examples
#' \dontrun{
#' lapply(seq(0, 1, 0.01), roc, vals, group = c("sexrate", "sample"), count.na = TRUE) %>%
#'   bind_rows %>%
#'   ungroup() %>%
#'   select(-sexrate) %>%
#'   spread(score, positive_fraction) -> y
#' ggplot(y, aes(x = Miss, y = Hit, color = sample)) +
#'   geom_line() +
#'   geom_point() +
#'   geom_abline(slope = 1, lty = 2)
#' }
roc <- function(alpha = 0.05, df, compare = c("0.0000", "1.0000"),
                stat = "p.rD", augment = NULL, count.na = TRUE,
                group = c("sexrate", "run", "seed", "sample")){

  pf <- lazyeval::interp(~n_hits(stat, alpha, count.na)/n(),
                         stat = as.name(stat), alpha = alpha, count.na = count.na)

  # pfsd <- list(positive_fraction_sd = ~sqrt((positive_fraction * (1 - positive_fraction))/ (n() - 1)))
  score_column <- list(score = ~ifelse(sexrate == compare[1], "True Positive", "False Positive"))
  df %>%
    filter_(.dots = list(~sexrate %in% compare)) %>%
    mutate_(.dots = list(ntp = ~sum(sexrate == compare[1]),
                         nfp = ~sum(sexrate == compare[2]))) %>%
    group_by_(.dots = c(group, "ntp", "nfp")) %>%
    summarize_(.dots = list(positive_fraction = pf, alpha = alpha)) %>%
    mutate_(.dots = score_column) %>%
    ungroup() %>%
    mutate_(.dots = list(sexrate = ~compare[1]))
}

n_hits <- function(stat, alpha, count.na = TRUE){
  res <- sum(stat <= alpha, na.rm = TRUE)
  res <- if (count.na) res + sum(is.na(stat)) else res
  return(res)
}

