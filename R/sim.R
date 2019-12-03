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

#' @title Simulate data with naturalistic autocorrelated errors 
#'
#' @param n_subj Number of subjects to simulate.
#' @param fixed Vector of population parameters for the fixed effects (intercept, main effect of A, main effect of B, AB interaction).
#' @param re_range Two-element vector defining the range for random effect variance (intercept and slope of A).
#' @details Simulates data from a 2x2 mixed design, with factor A
#'   within subjects and factor B between subjects, with residuals
#'   randomly sampled from the KKL dataset (see ?KKL in the
#'   RePsychLing package).
#' @seealso \code{\link{kkl_df}}, \code{\link{kkl_mx}}
#' @return A data frame, with:
#' \describe{
#'    \item{ts_r}{the trial number given random ordering}
#'    \item{ts_b}{the trial number given blocked allocation of trials by the within-factor "A"}
#'    \item{Y_fit}{the 'fitted' value (combining fixed and random effects but not residual error)}
#'    \item{Y_acr}{the response variable (random order)}
#'    \item{Y_acb}{the response variable (blocked order)}
#' }
#' @importFrom magrittr %>%
#' @export
sim_2x2_kkl <- function(n_subj, fixed, re_range) {
  n_per <- 800L # number of trials

  design_args <- list(ivs = c(A = 2, B = 2),
                      n_item = n_per * 2L, 
                      between_item = c("A", "B"),
                      between_subj = c("B"))

  vr <- if (length(re_range) == 1L) rep(re_range, 2L) else re_range
  
  parms <- funfact::gen_pop(design_args, n_subj, var_range = vr)
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- 0

  if (n_subj > nrow(kkl_mx))
    stop("n_subj cannot exceed number of KKL subjects (86)")

  if (n_subj %% 2)
    stop("n_subj must be divisible by 2")
      
  dat <- funfact::sim_norm(design_args, n_subj, parms, verbose = TRUE) %>%
    dplyr::mutate(Y_fit = Y - err)

  res <- c(t(kkl_mx[sample(nrow(kkl_mx), n_subj), ]))
  
  ## trials in random order
  dat2 <- dat %>%
    dplyr::group_by(subj_id) %>%
    dplyr::mutate(ts_r = sample(seq_len(n_per))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_r) %>%
    dplyr::mutate(Y_acr = Y_fit + res)

  ## trials blocked by level of A (randomized within)
  dat2 %>%
    dplyr::group_by(subj_id, A) %>%
    dplyr::mutate(ts_b = sample(seq_len(n_per / 2) +
                                (A == "A2") * (n_per / 2))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_b) %>%
    dplyr::mutate(Y_acb = Y_fit + res) %>%
    dplyr::select(subj_id, item_id, ts_r, ts_b, A, B,
                  Y_fit, Y_acr, Y_acb)
}
