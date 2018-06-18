#' Plot autocorrelation function
#'
#' Estimate and plot the autocorrelation function.
#'
#' @inheritParams acf_subj
#' @importFrom rlang "!!"
#' @export
acf_plot <- function(data, subject, residuals, lag.max = NULL) {
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  acorr <- acf_subj(data, !!subj, !!resid, lag.max)
  acorr_ci <- acf_ci(data, !!subj, !!resid, lag.max)
  if (is.null(lag.max)) {
    lag.max <- max(acorr_ci[["lag"]])
  }
  ggplot2::ggplot(acorr, ggplot2::aes(lag)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), acorr_ci,
                         fill = "red", alpha = .2) +
      ggplot2::geom_line(ggplot2::aes(y = r, group = subj), alpha = .2)
}
