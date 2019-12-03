#' Simulate data with sine-wave autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate.
#' @param n_obs Number of observations per subject.
#' @param fixed Vector of population parameters for the fixed effects (intercept, main effect of A, main effect of B, AB interaction).
#' @param err_range Two-element vector defining the range for sigma^2, the residual error parameter.
#' @param re_range Two-element vector defining the range for random effect variance (intercept and slope of A).
#' @param phase Whether the phase of each participant's time-varying function is offset by a random angle (from -pi to pi).
#' @param amp Whether the amplitude of each participant's time-varying function is modulated by a random value (from 0 to 2).
#' @param verbose Whether to include random intercept, random slope, phase angle, and amplitude for each participant.
#' @details Simulates data from a 2x2 mixed design, with factor A
#'   within subjects and factor B between subjects.
#' @return A data frame, with \code{ts_r} as the trial number ordered
#'   randomly; \code{ts_b} as the trial number blocked by the within
#'   factor ("A"); \code{Y_fit} as the 'fitted' value (combining fixed
#'   and random effects but not residual error); \code{Y_acr} as the
#'   response variable with autocorrelation; and \code{Y_acb} as the
#'   response variable with autocorrelation for the blocked data.
#' @importFrom magrittr %>%
#' @export
sim_2x2_sin <- function(n_subj, n_obs, fixed,
                        err_range, re_range,
                        phase = TRUE, amp = FALSE,
                        verbose = FALSE) {
  x <- seq(-pi, pi, length.out = n_obs)

  step_begin <- if (phase) runif(n_subj, -pi, pi) else rep(0, n_subj)
  step_amp <- if (amp) runif(n_subj, 0, 2) else rep(1, n_subj)

  err <- runif(1, err_range[1], err_range[2])
  
  resids <- purrr::map2(step_begin, step_amp,
                        ~ rnorm(n_obs, .y * (sin(x + .x) / sd(sin(x))),
                                sqrt(err)))

  rtbl <- tibble::tibble(subj_id = factor(seq_len(n_subj)), phase = step_begin,
                         amp = step_amp) %>%
    dplyr::mutate(
             resid = purrr::map2(
                              step_begin, step_amp,
                              ~ tibble::tibble(t = seq_along(x),
                                               x = x,
                                               sin = .y * (sin(x + .x) / sd(sin(x))),
                                               err = rnorm(n_obs, 0, sqrt(err)),
                                               resid = sin + err))) %>%
    tidyr::unnest(c(resid))
  
  my_design <- list(ivs = c(A = 2, B = 2),
		    n_item = n_obs * 2L,
		    between_item = c("A", "B"),
		    between_subj = c("B"))

  parms <- funfact::gen_pop(my_design, n_subj, var_range = re_range)
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- 0

  dat <- funfact::sim_norm(my_design, n_subj, parms, verbose = TRUE) %>%
    dplyr::mutate(subj_id = factor(subj_id),
		  list_id = factor(list_id),
		  item_id = factor(item_id)) %>%
    dplyr::rename(Y_fit = Y) %>%
    dplyr::select(-err, -ire)

  n_per <- dat %>% dplyr::count(subj_id) %>% dplyr::pull(n) %>% unique()
  stopifnot(length(n_per) == 1L)
  
  ## trials in random order
  dat2 <- dat %>%
    dplyr::group_by(subj_id) %>%
    dplyr::mutate(ts_r = sample(seq_len(n_per))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_r) %>%
    dplyr::inner_join(rtbl, c("subj_id" = "subj_id",
                              "ts_r" = "t")) %>%
    dplyr::rename(x_r = x, sin_r = sin, err_r = err, res_r = resid) %>%
    dplyr::mutate(Y_acr = Y_fit + res_r)

  ## trials blocked by level of A (randomized within)
  dat3 <- dat2 %>%
    dplyr::group_by(subj_id, A) %>%
    dplyr::mutate(ts_b = sample(seq_len(n_per / 2) +
				(A == "A2") * (n_per / 2))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(subj_id, ts_b) %>%
    dplyr::inner_join(rtbl %>% dplyr::select(-phase, -amp),
                      c("subj_id" = "subj_id",
                        "ts_b" = "t")) %>%
    dplyr::rename(x_b = x, sin_b = sin, err_b = err, res_b = resid) %>%  
    dplyr::mutate(Y_acb = Y_fit + res_b,
		  tnum_r = (ts_r - (n_obs + 1) / 2) / (n_obs - 1),
		  tnum_b = (ts_b - (n_obs + 1) / 2) / (n_obs - 1)) %>%
    funfact::with_dev_pred(c("A", "B")) %>%
    dplyr::select(subj_id, A, B, AA2, BB2, tnum_b, tnum_r, Y_acb, Y_acr,
           ts_b, x_b, sin_b, err_b,
           ts_r, x_r, sin_r, err_r,
           Y_fit, fix_y, sre, phase, amp)

  if (verbose) {
    ## get random effects
    rfx <- dplyr::distinct(dat, subj_id, A, sre) %>%
      tidyr::spread(A, sre) %>%
      dplyr::mutate(sri = (A1 + A2) / 2,#
                    srs = (A2 - A1)) %>%
      dplyr::select(-A1, -A2)

    dat3 %>%
      dplyr::inner_join(rfx, "subj_id") %>%
      dplyr::select(subj_id:fix_y, phase, amp, sre, sri, srs)
  } else {
    dat3 %>%
      dplyr::select(subj_id:Y_acr)
  }
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
  my_formula_rhs_b <- my_formula_rhs # blocked version
  
  my_formula_rhs_no_gamm <- paste0(my_formula_rhs,
                                   " + \n   s(subj_id, bs = \"re\")")
  
  if (cs) {
    my_formula_rhs <- paste0(my_formula_rhs,
                             " + \n   s(tnum_r, bs = \"tp\")")
    my_formula_rhs_b <- paste0(my_formula_rhs_b,
                             " + \n   s(tnum_b, bs = \"tp\")")
  }

  if (by_subj_fs) {
    my_formula_rhs <- paste(my_formula_rhs,
                            " + \n   s(tnum_r, subj_id, bs = \"fs\")")
    my_formula_rhs_b <- paste(my_formula_rhs_b,
                              " + \n   s(tnum_b, subj_id, bs = \"fs\")")
  }

  ## fit the GAMM models using mgcv::bam
  mod_rand <- mgcv::bam(as.formula(paste0("Y_acr ~", my_formula_rhs)),
                        data = dat)
  mod_block <- mgcv::bam(as.formula(paste0("Y_acb ~", my_formula_rhs_b)),
                         data = dat)

  ## fit the LMEM models using mgcv::bam
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
}

#' Monte Carlo Simulation of Sine Wave Autocorrelation
#'
#' @param nmc Number of Monte Carlo runs.
#' @param nsubj Number of subjects per dataset.
#' @param ntrials Number of trials per subject.
#' @param fixed Fixed effect size (four-element vector with intercept, main effect of A, main effect of B, and AB interaction).
#' @param err_range Error variance (sigma^2).
#' @param re_range Range of random effects variance.
#' @param varying_phase Does the sine-wave autocorrelation have varying (vs fixed) phase over subjects?
#' @param varying_amp Does the sine-wave autocorrelation have varying (vs fixed) amplitude over subjects?
#' @return An array with labeled dimensions of class 'sinemcs'.
#' @export
sine_sim <- function(nmc,
                     nsubj = 48L,
                     ntrials = 48L,
                     fixed = rep(0, 4),
                     err_range = c(1, 3),
                     re_range = c(1, 3),
                     varying_phase = TRUE,
                     varying_amp = FALSE) {
  rmx <- replicate(nmc, {
    dat <- sim_2x2_sin(nsubj, ntrials,
                       fixed,
                       err_range = err_range,
                       re_range = re_range,
                       phase = varying_phase,
                       amp = varying_amp)
    fit_2x2_sin(dat, cs = !varying_phase)
  })
  class(rmx) <- "sinemcs"
  rmx
}

#' @export
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

#' @export
print.sinemcs <- function(x) {
  summary(x)
}
