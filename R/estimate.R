## functions in this file estimate autocorrelation functions

#' @title Estimate by-subject acfs
#'
#' Estimate autocorrelation functions for each subject in a dataset.
#'
#' @param subject Unquoted variable name in \code{data} identifying
#'   individual subjects.
#' @param residuals Unquoted variable name in \code{data} of
#'   residuals.
#' @param lag.max maximum lag at which to calculate the acf.  Default
#'   is to run two passes over the dataset. In the first pass, the
#'   length of each subject's autocorrelation function is given by
#'   10*log(N, 10), where N is the number of observations in the
#'   series.  In the second pass, the length of the longest series
#'   will be identified, and this will be used for \code{lag.max}.
#' @return A table with the autocorrelation function for each subject.
#' @seealso \code{\link{acf}}
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @export
acf_subj <- function(data, subject, residuals, lag.max = NULL) {
  acf_tbl <- function(x, col, lag.max) {
    vec <- dplyr::pull(x, !!col)
    res <- c(acf(vec, lag.max = lag.max, plot = FALSE)[[1]])
    tibble::tibble(lag = seq_along(res), r = res)
  }
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  dnest <- data %>%
    dplyr::select(!!subj, !!resid) %>%
    dplyr::group_by(!!subj) %>%
    tidyr::nest(.key = "res")
  acf_df <- dnest %>%
    dplyr::mutate(acf = purrr::map(res, acf_tbl, resid, lag.max)) %>%
    dplyr::select(subject, acf) %>%
    tidyr::unnest()
  if (is.null(lag.max)) {
    max_lag <- acf_df %>%
      dplyr::group_by(subject) %>%
      dplyr::summarize(max_lag = max(lag)) %>%
      dplyr::pull(max_lag) %>% max()
    message("setting lag.max to ", max_lag)
    acf_df <- dnest %>%
      dplyr::mutate(acf = purrr::map(res, acf_tbl, resid, max_lag)) %>%
      dplyr::select(subject, acf) %>%
      tidyr::unnest()
  }
  acf_df
}

#' @title Confidence interval for autocorrelation function
#'
#' Calculate a permutation-based confidence interval for the
#' by-subject autocorrelation function
#'
#' @param data A dataset.
#' @param subject Unquoted name of the variable identifying individual
#'   subjects.
#' @param residuals Unquoted variable name for residuals.
#' @param lag.max Maximum lag at which to calculate the autocorrelation function. See \code{\link{acf_subj}} for default behavior.
#' @param alpha Alpha level.
#' @details Permutes the residuals for each participant before calculating the confidence interval for the correlation at each time lag.
#' @return Table with (1-alpha)% confidence interval for the autocorrelation function at each time lag.
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @export
acf_ci <- function(data, subject, residuals, lag.max = NULL, alpha = .05) {
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  if (is.null(lag.max)) {
    lag.max <- acf_subj(data, !!subj, !!resid) %>%
      dplyr::pull(lag) %>% max()
  }
  permute_resids(data, !!subj, !!resid) %>%
    acf_subj(!!subj, !!resid, lag.max) %>%
    dplyr::group_by(lag) %>%
    dplyr::summarize(lower = quantile(r, probs = alpha / 2),
		     upper = quantile(r, probs = 1 - (alpha / 2)))
}

#' @title Durbin-Watson statistic by subject
#'
#' Calculates the Durbin-Watson statistic by subject.
#'
#' @param data A dataset.
#' @param subject Unquoted variable name in \code{data} identifying
#'   individual subjects.
#' @param residuals Unquoted variable name in \code{data} of
#'   residuals.
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#' @export
dw_subj <- function(data, subject, residual) {
  dw_calc <- function(x) {
    sum((x[-1] - x[-length(x)])^2) / sum(x^2)
  }
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residual)
  data %>%
    dplyr::group_by(!!subj) %>%
    dplyr::summarize(dw = dw_calc(!!resid)) %>%
    dplyr::ungroup()
}
