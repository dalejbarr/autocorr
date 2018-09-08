#' Simulate data with sine-wave autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate.
#' @param n_obs Number of observations per subject.
#' @param fixed Vector of population parameters for the fixed effects (intercept, main effect of A, main effect of B, AB interaction).
#' @param err Sigma, the residual error parameter.
#' @param phase Whether the phase of each participant's time-varying function is offset by a random angle (from -pi to pi).
#' @param amp Whether the amplitude of each participant's time-varying function is modulated by a random value (from 0 to 2).
#' @param re_varmax Maximum variance for random effects parameters (slopes and intercepts).
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
sim_2x2_sin <- function(n_subj, n_obs, fixed, err,
                        phase = TRUE, amp = FALSE,
                        re_varmax = 3) {
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

  parms <- gen_pop(my_design, n_subj, var_range = c(0, re_varmax))
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- err

  dat <- sim_norm(my_design, n_subj, parms, verbose = TRUE) %>%
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
    with_dev_pred(c("A", "B"))
}

#' Fit models to 2x2 data with sine wave errors
#'
#' @param dat Data generated using \link{sim_2x2_sin}.
#' @param cs Fit gamm smooth as main effect of trial (in addition to by-subject factor smooths).
#' @param by_subj_fs Include by-subject factor smooths?
#' @return Array with parameter estimates and standard errors.
#' @export
fit_2x2_sin <- function(dat, cs = FALSE, by_subj_fs = TRUE) {
  ## function to extract model statistics
  mod_stats <- function(m_gamm, m_lmem) {
    mc <- anova(m_gamm, m_lmem)
    m_gamm_s <- summary(m_gamm)
    m_lmem_s <- summary(m_lmem)
    v <- c(coef(m_gamm)[1:4],
           sqrt(diag(vcov(m_gamm)[1:4, 1:4])),
           m_gamm_s[["p.table"]][, "Pr(>|t|)"],
           AIC(m_gamm), mc[["Resid. Df"]][1],
           mc[["Resid. Dev"]][1],
           coef(m_lmem)[1:4],
           sqrt(diag(vcov(m_lmem)[1:4, 1:4])),
           m_lmem_s[["p.table"]][, "Pr(>|t|)"],
           AIC(m_lmem), mc[["Resid. Df"]][2],
           mc[["Resid. Dev"]][2])
    vn <- c("e_int", "e_A", "e_B", "e_AB",
            "se_int", "se_A", "se_B", "se_AB",
            "p_int", "p_A", "p_B", "p_AB",
            "AIC", "resid_df", "resid_dev")
    names(v) <- rep(vn, 2)
    v
  }
  ## determine the model formula
  my_formula_rhs <- "AA2 * BB2 + \n   s(subj_id, AA2, bs = \"re\")"
  my_formula_rhs_no_gamm <- paste0(my_formula_rhs,
                                   " + \n   s(subj_id, bs = \"re\")")
  
  if (cs) {
    my_formula_rhs <- paste0(my_formula_rhs,
                             " + \n   s(tnum_r, bs = \"tp\")")
  }

  if (by_subj_fs) {
    my_formula_rhs <- paste(my_formula_rhs,
                            " + \n   s(tnum_r, subj_id, bs = \"fs\")")
  }

  ## fit the GAMM models using mgcv::bam
  mod_rand <- mgcv::bam(as.formula(paste0("Y_acr ~", my_formula_rhs)),
                        data = dat)
  mod_block <- mgcv::bam(as.formula(paste0("Y_acb ~", my_formula_rhs)),
                         data = dat)

  ## fit the non-GAMM models using mgcv::bam
  mod_rand_2 <- mgcv::bam(as.formula(paste0("Y_acr ~", my_formula_rhs_no_gamm)),
                        data = dat)
  mod_block_2 <- mgcv::bam(as.formula(paste0("Y_acb ~", my_formula_rhs_no_gamm)),
                        data = dat)

  array(c(mod_stats(mod_rand, mod_rand_2),
          mod_stats(mod_block, mod_block_2)),
        dim = c(15, 2, 2),
        dimnames = list(parm = names(mod_stats(mod_rand, mod_rand_2))[1:15],
                        mod = c("GAMM", "LMEM"),
                        vers = c("randomized", "blocked")))  
  ## array(c(coef(mod_no)[1:4],
  ##  sqrt(diag(vcov(mod_no)[1:4, 1:4])),
  ##  coef(mod_rand)[1:4],
  ##  sqrt(diag(vcov(mod_rand)[1:4, 1:4])),
  ##  coef(mod_block)[1:4],
  ##  sqrt(diag(vcov(mod_block)[1:4, 1:4])),
  ##  coef(mod_no_2)[1:4],
  ##  sqrt(diag(vcov(mod_no_2)[1:4, 1:4])),
  ##  coef(mod_rand_2)[1:4],
  ##  sqrt(diag(vcov(mod_rand_2)[1:4, 1:4])),
  ##  coef(mod_block_2)[1:4],
  ##  sqrt(diag(vcov(mod_block_2)[1:4, 1:4]))),
  ##  dim = c(4, 2, 6),
  ##  dimnames = list(coef = c("(Intercept)", "A", "B", "A:B"),
  ## 	param = c("est", "stderr"),
  ## 	model = c("gamm_no", "gamm_rand", "gamm_block",
  ## 		  "nogamm_no", "nogamm_rand", "nogamm_block")))
}

#' Monte Carlo Simulation of Sine Wave Autocorrelation
#'
#' @param nmc Number of Monte Carlo runs.
#' @param nsubj Number of subjects per dataset.
#' @param ntrials Number of trials per subject.
#' @param fixed Fixed effect size (four-element vector with intercept, main effect of A, main effect of B, and AB interaction).
#' @param err Error variance (sigma^2).
#' @param varying_phase Does the sine-wave autocorrelation have varying (vs fixed) phase over subjects?
#' @param varying_amp Does the sine-wave autocorrelation have varying (vs fixed) amplitude over subjects?
#' @param re_varmax Maximum of the range for the random effect variance parameters.
#' @return An array with labeled dimensions of class 'sinemcs'.
#' @export
sine_sim <- function(nmc,
                     nsubj = 48L,
                     ntrials = 48L,
                     fixed = rep(0, 4),
                     err = 6,
                     varying_phase = TRUE,
                     varying_amp = FALSE,
                     re_varmax = 3) {
  rmx <- replicate(nmc,
                   fit_2x2_sin(
                     sim_2x2_sin(nsubj, ntrials,
                                 fixed, err,
                                 phase = varying_phase,
                                 amp = varying_amp,
                                 re_varmax = re_varmax),
                     cs = !varying_phase))
  class(rmx) <- "sinemcs"
  rmx
}

summary.sinemcs <- function(x, alpha = .05) {
  cat("Sine Wave Autocorrelation Simulation Results\n\n")
  cat("Number of Monte Carlo Runs: ", dim(x)[4], "\n\n")
  cat("Proportion Significant:\n")
  pvals <- x[c("p_A", "p_B", "p_AB"), , , ]
  sig <- pvals < alpha
  nsig <- apply(sig, 1:3, sum, na.rm = TRUE)
  ntot <- apply(sig, 1:3, function(.x) sum(!is.na(.x)))
  cat("Randomized:\n")
  print((nsig / ntot)[, , 1])
  cat("\nBlocked:\n")
  print((nsig / ntot)[, , 2])
}

print.sinemcs <- function(x) {
  summary(x)
}
