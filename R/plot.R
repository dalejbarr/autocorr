#' Plot autocorrelation function
#'
#' Estimate and plot the autocorrelation function.
#'
#' @inheritParams acf_subj
#' @param log_log Whether to produce a log-log plot.
#' @importFrom rlang "!!"
#' @export
acf_plot <- function(data, subject, residuals, lag.max = NULL,
                     log_log = FALSE) {
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  acorr <- acf_subj(data, !!subj, !!resid, lag.max)
  acorr_ci <- acf_ci(data, !!subj, !!resid, lag.max)
  if (is.null(lag.max)) {
    lag.max <- max(acorr_ci[["lag"]])
  }
  if (log_log) {
    ggplot2::ggplot(acorr, ggplot2::aes(log(lag))) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = log(lower), ymax = log(upper)),
                           acorr_ci,
                           fill = "red", alpha = .2) +
        ggplot2::geom_line(ggplot2::aes_(y = ~log(r), group = subj),
                           alpha = .2)    
  } else {
    ggplot2::ggplot(acorr, ggplot2::aes(lag)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), acorr_ci,
                           fill = "red", alpha = .2) +
        ggplot2::geom_line(ggplot2::aes_(y = ~r, group = subj),
                           alpha = .2)
  }
}

#' Plot Durbin-Watson statistics
#'
#' Calculate Durbin-Watson statistics for original (and permuted) residuals in the dataset and produce a density plot.
#'
#' @inheritParams acf_subj
#' @importFrom rlang "!!"
#' @importFrom magrittr "%>%"
#' @export
dw_plot <- function(data, subject, residuals, lag.max = NULL) {
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  dw_o <- dw_subj(data, !!subj, !!resid)
  dw_p <- dw_subj(permute_resids(data, !!subj, !!resid), !!subj, !!resid)

  dplyr::bind_rows(dw_o %>% dplyr::mutate(residual = "original"),
                   dw_p %>% dplyr::mutate(residual = "permuted")) %>%
    ggplot2::ggplot(ggplot2::aes(dw, fill = residual)) +
      ggplot2::geom_density(alpha = .2) +
      ggplot2::coord_cartesian(xlim = c(0, 4))
}
