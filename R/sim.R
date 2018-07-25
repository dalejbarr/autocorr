#' @title Simulate autocorrelated residual error based on source function
#'
#' @param source_acf Source autocorrelation function to use.
#' @param length.out Length of the output time series.
#' @param sigma Desired standard deviation of output time series.
#' @param center Mean center the time series.
#' @return A vector of residuals modeled after \code{source_acf}.
#' @export
sim_acerr <- function(source_acf, length.out = 48L, sd = 1, 
                      center = TRUE) {
  num_nas <- length(source_acf) - 1L
  nas_head <- floor(num_nas / 2) + 1L
  vec <- stats::filter(rnorm(length.out + num_nas, sd = sd), rev(source_acf))
  res <- vec[seq(nas_head, by = 1L, length.out = length.out)]
  if (center)
    res - mean(res)
  else
    res
}

#' @title Simulate data with autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate.
#' @param n_obs Number of observations per subject.
#' @param params Population parameter values for the data generating
#'   process; the structure should be identical to that returned by
#'   \code{\link{funfact::gen_pop}}.
#' @param is_acf Is argument \code{amr} a list of residuals (TRUE) or a matrix
#'   of autocorrelation functions (FALSE)?
#' @param amx List of residuals or matrix of autocorrelation functions, with each element or row representing a different subject.
#' @param amx_wt Weights for sampling elements/rows from \code{amx}; NULL gives equal probability.
#' @details Simulates data from a 2x2 mixed design, with factor A
#'   within subjects and factor B between subjects.
#' @return A data frame, with \code{ts_r} as the trial number ordered
#'   randomly; \code{ts_b} as the trial number blocked by the within
#'   factor ("A"); \code{Y_fit} as the 'fitted' value (combining fixed
#'   and random effects but not residual error); \code{Y_acn} as the
#'   response variable without autocorrelation; \code{Y_acr} as the
#'   response variable with autocorrelation; and \code{Y_acb} as the
#'   response variable with autocorrelation for the blocked data.
#' @importFrom magrittr %>%
#' @export
sim_2x2 <- function(n_subj, n_obs, params, is_acf,
                    amx, amx_wt = NULL) {
  design_args <- list(ivs = c(A = 2, B = 2),
                      n_item = n_obs * 2L,
                      between_item = c("A", "B"),
                      between_subj = c("B"))
  dat <- funfact::sim_norm(design_args, n_subj, params, verbose = TRUE) %>%
    dplyr::mutate(Y_fit = Y - err)

  n_per <- dat %>% dplyr::count(subj_id) %>% dplyr::pull(n) %>% unique()
  stopifnot(length(n_per) == 1L)

  if (is_acf) {
    if (!is.matrix(amx)) {
      stop("is_acf was TRUE but amx was not a matrix")
    }
    acerr_all <- purrr::map(sample(seq_len(nrow(amx)), n_subj, TRUE,
                                   prob = amx_wt),
                            ~ sim_acerr(amx[.x, ], n_per, params$err_var))
  } else {
    if (!is.list(amx)) {
      stop("is_acf was FALSE but amx was not a list")
    }
    if (any(purrr::map_int(amx, length) != n_per)) {
      stop("elements of 'amx' must be vectors of same length as number of observations for a given subject (", n_per, ")")
    }
    acerr_all <- amx[sample(seq_len(length(amx)), n_subj, TRUE, prob = amx_wt)]
  }

  acerr_no_ac <- purrr::map(acerr_all, sample)

  ## trials in random order
  dat2 <- dat %>%
    dplyr::group_by(subj_id) %>%
    dplyr::mutate(ts_r = sample(seq_len(n_per))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_r) %>%
    dplyr::mutate(Y_acr = Y_fit + unlist(acerr_all),
                  Y_acn = Y_fit + unlist(acerr_no_ac))

  ## trials blocked by level of A (randomized within)
  dat2 %>%
    dplyr::group_by(subj_id, A) %>%
    dplyr::mutate(ts_b = sample(seq_len(n_per / 2) +
                                (A == "A2") * (n_per / 2))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_b) %>%
    dplyr::mutate(Y_acb = Y_fit + unlist(acerr_all)) %>%
    dplyr::select(subj_id, item_id, ts_r, ts_b, A, B,
                  Y_fit, Y_acn, Y_acr, Y_acb)
}
