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
roc <- function(alpha = 0.05, df, compare = c("0.0000", "1.0000"), stat = "p.rD",
                count.na = TRUE, group = c("sexrate", "run", "seed", "sample")){
  if (count.na){
      pf <- lazyeval::interp(~(sum(stat <= alpha, na.rm = TRUE) + sum(is.na(stat)))/n(),
                         stat = as.name(stat), alpha = alpha)
  } else {
      pf <- lazyeval::interp(~sum(stat <= alpha, na.rm = TRUE)/n(),
                         stat = as.name(stat), alpha = alpha)
  }
  transform <- lazyeval::interp(~ifelse(sexrate == ~x, "Hit", "Miss"),
                                x = compare[1])
  df %>%
    filter_(.dots = list(~sexrate %in% compare)) %>%
    group_by_(.dots = group) %>%
    summarize_(.dots = list(positive_fraction = pf,
                            alpha = alpha)) %>%
    mutate_(.dots = list(score = ~ifelse(sexrate == compare[1], "True Positive", "False Positive"))) %>%
    ungroup() %>%
    mutate_(.dots = list(sexrate = compare[1]))
}

roc_curve <- function(alpha = seq(0, 1, by = 0.05), df, ...){

}
