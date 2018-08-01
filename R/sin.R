#' Simulate data with sine-wave autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate.
#' @param n_obs Number of observations per subject.
#' @param fixed Vector of population parameters for the fixed effects (intercept, main effect of A, main effect of B, AB interaction).
#' @param err Sigma, the residual error parameter.
#' @param phase Whether the phase of each participant's time-varying function is offset by a random angle (from -pi to pi).
#' @param amp Whether the amplitude of each participant's time-varying function is modulated by a random value (from 0 to 2).
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
sim_2x2_sin <- function(n_subj, n_obs, fixed, err, phase = TRUE, amp = FALSE) {
  x <- seq(-pi, pi, length.out = n_obs)

  step_begin <- if (phase) runif(n_subj, -pi, pi) else rep(0, n_subj)
  step_amp <- if (amp) runif(n_subj, 0, 2) else rep(1, n_subj)

  resids <- purrr::map2(step_begin, step_amp,
                        ~ rnorm(n_obs, .y * (sin(x + .x) / sd(sin(x))), err))
  resids_no_ar <- purrr::map(resids, sample) # shuffled

  my_design <- list(ivs = c(A = 2, B = 2),
		    n_item = n_obs * 2L,
		    between_item = c("A", "B"),
		    between_subj = c("B"))

  parms <- funfact::gen_pop(my_design, n_subj)
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- err

  dat <- funfact::sim_norm(my_design, n_subj, parms, verbose = TRUE) %>%
    dplyr::mutate(subj_id = factor(subj_id),
		  list_id = factor(list_id),
		  item_id = factor(item_id),
		  Y_fit = Y - err)

  n_per <- dat %>% dplyr::count(subj_id) %>% dplyr::pull(n) %>% unique()
  stopifnot(length(n_per) == 1L)

  ## trials in random order
  dat2 <- dat %>%
    dplyr::group_by(subj_id) %>%
    dplyr::mutate(ts_r = sample(seq_len(n_per))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_r) %>%
    dplyr::mutate(Y_acr = Y_fit + unlist(resids),
		  Y_acn = Y_fit + unlist(resids_no_ar))

  ## trials blocked by level of A (randomized within)
  dat2 %>%
    dplyr::group_by(subj_id, A) %>%
    dplyr::mutate(ts_b = sample(seq_len(n_per / 2) +
				(A == "A2") * (n_per / 2))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_b) %>%
    dplyr::mutate(Y_acb = Y_fit + unlist(resids),
		  tnum_r = (ts_r - (n_obs + 1) / 2) / (n_obs - 1),
		  tnum_b = (ts_b - (n_obs + 1) / 2) / (n_obs - 1)) %>%
    dplyr::select(subj_id, item_id, ts_r, ts_b, A, B,
		  Y_fit, Y_acn, Y_acr, Y_acb, tnum_r, tnum_b) %>%
    funfact::with_dev_pred(c("A", "B"))
}

#' Fit models to 2x2 data with sine wave errors
#'
#' @param dat Data generated using \link{sim_2x2_sin}.
#' @return Array with parameter estimates and standard errors.
#' @export
fit_2x2_sin <- function(dat) {
  ## fit the GAMM models using mgcv::bam
  mod_no <- mgcv::bam(Y_acn ~ AA2 * BB2 +
			s(tnum_r, subj_id, bs = "fs") +
			s(subj_id, AA2, bs = "re"),
		      data = dat)

  mod_rand <- mgcv::bam(Y_acr ~ AA2 * BB2 +
			  s(tnum_r, subj_id, bs = "fs") +
			  s(subj_id, AA2, bs = "re"),
			data = dat)

  mod_block <- mgcv::bam(Y_acb ~ AA2 * BB2 +
			   s(tnum_b, subj_id, bs = "fs") +
			   s(subj_id, AA2, bs = "re"),
			 data = dat)

  ## fit the non-GAMM models using mgcv::bam
  mod_no_2 <- mgcv::bam(Y_acn ~ AA2 * BB2 +
			  s(subj_id, bs = "re") +
			  s(subj_id, AA2, bs = "re"),
			data = dat)

  mod_rand_2 <- mgcv::bam(Y_acr ~ AA2 * BB2 +
			    s(subj_id, bs = "re") +
			    s(subj_id, AA2, bs = "re"),
			  data = dat)

  mod_block_2 <- mgcv::bam(Y_acb ~ AA2 * BB2 +
			     s(subj_id, bs = "re") +
			     s(subj_id, AA2, bs = "re"),
			   data = dat)

  array(c(coef(mod_no)[1:4],
	  sqrt(diag(vcov(mod_no)[1:4, 1:4])),
	  coef(mod_rand)[1:4],
	  sqrt(diag(vcov(mod_rand)[1:4, 1:4])),
	  coef(mod_block)[1:4],
	  sqrt(diag(vcov(mod_block)[1:4, 1:4])),
	  coef(mod_no_2)[1:4],
	  sqrt(diag(vcov(mod_no_2)[1:4, 1:4])),
	  coef(mod_rand_2)[1:4],
	  sqrt(diag(vcov(mod_rand_2)[1:4, 1:4])),
	  coef(mod_block_2)[1:4],
	  sqrt(diag(vcov(mod_block_2)[1:4, 1:4]))),
	dim = c(4, 2, 6),
	dimnames = list(coef = c("(Intercept)", "A", "B", "A:B"),
			param = c("est", "stderr"),
			model = c("gamm_no", "gamm_rand", "gamm_block",
				  "nogamm_no", "nogamm_rand", "nogamm_block")))
}
