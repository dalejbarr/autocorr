#' Fit logistic regression model using lmer
#'
#' Convenience function that fits a logit model and returns NULL when the estimation process fails to converge.
#' @param formula Model formula.
#' @param data Data frame containing variables referenced in formula.
#' @return If estimation process converges, a \code{glmerMod} object; otherwise \code{NULL}.
#' @export
fit_logit_lmer <- function(formula, data) {
  tryCatch(glmer(formula, data, binomial),
           error = function(e) NULL,
           warning = function(w) NULL)
}

#' Fit model via Julia MixedModels package
#'
#' Fits a linear mixed-effects model via the Julia MixedModels
#' package.
#' 
#' @param mform Model formula.
#' @param dat Data frame.
#' @param jsess Julia session (established by \code{JuliaCall::julia_setup()}).
#' @param name Name of the Julia object to create.
#' @return A Julia fitted model object.
#' @export
fit_jl <- function(formula, data, jsess, name = "fm1") {
  jsess$assign("form", formula)
  jsess$assign("dat", data)

  jsess$eval(paste0(name, " = fit(LinearMixedModel, form, dat)"))
  jsess$eval(name)
}

#' Fit logit model via Julia MixedModels package
#'
#' Fits a mixed-effects logistic regression model via the Julia MixedModels
#' package.
#' 
#' @param mform Model formula.
#' @param dat Data frame.
#' @param jsess Julia session (established by \code{JuliaCall::julia_setup()}).
#' @param name Name of the Julia object to create.
#' @return A Julia fitted model object.
#' @export
fit_logit_jl <- function(formula, data, jsess, name = "fm1") {
  jsess$assign("form", formula)
  jsess$assign("dat", data)
  jsess$eval(paste0(name, " = fit(GeneralizedLinearMixedModel, form, dat, Bernoulli())"))
  jsess$eval(name)
}


#' Simulate and fit 2x2 autocorrelated data
#'
#' Simulate data from a 2x2 mixed design (A within and B between) with
#' autocorrelated residuals.
#'
#' @param n_subj Number of subjects.
#' @param n_obs Number of observations per subject.
#' @param fixed vector of fixed effects in the following order:
#'   Intercept, main effect of A, main effect of B, and AB
#'   interaction.
#' @param jsess Julia session.
#' @param amx Matrix of autocorrelation functions.
#' @param amx_wt Weights for sampling rows from \code{amx}.
#' @details Simulate data corresponding to three models: no
#'   autocorrelation (\code{mod_no}), random trial order
#'   (\code{mod_rand}), and with trials blocked by levels of A
#'   (\code{mod_block}). Models are fit via Julia (using
#'   \code{\link{run_jl}}).
#' @return A two dimensional array with fixed effects estimates and
#'   standard errors.
#' @seealso \code{\link{sim_2x2}}, \code{\link{fit_jl}}
#' @importFrom magrittr %>%
#' @export
run_2x2 <- function(n_subj, n_obs, fixed, jsess, amx, amx_wt = NULL) {
  dat <- gen_2x2(n_subj, n_obs, fixed, amx, amx_wt)

  ## fit the models using autocorr::fit_jl
  mod_no <- fit_jl(Y_acn ~ AA2 * BB2 + (AA2 | subj_id), dat,
                   jsess, "mod_no")
  mod_rand <- fit_jl(Y_acr ~ AA2 * BB2 + (AA2 | subj_id), dat,
                     jsess, "mod_rand")
  mod_block <- fit_jl(Y_acb ~ AA2 * BB2 + (AA2 | subj_id), dat,
                      jsess, "mod_block")

  array(c(jsess$eval("StatsBase.coef(mod_no)"),
          jsess$eval("StatsBase.stderror(mod_no)"),
          jsess$eval("StatsBase.coef(mod_rand)"),
          jsess$eval("StatsBase.stderror(mod_rand)"),
          jsess$eval("StatsBase.coef(mod_block)"),
          jsess$eval("StatsBase.stderror(mod_block)")),
        dim = c(4, 2, 3),
        dimnames = list(coef = c("(Intercept)", "A", "B", "A:B"),
                        param = c("est", "stderr"),
                        model = c("mod_no", "mod_rand", "mod_block")))
}

#' Simulate and Fit GAMM (OLD)
#'
#' Simulate 2x2 data and fit GAMM models (deprecated: use gamm_2x2)
#'
#' @param n_subj Number of subjects.
#' @param n_obs Number of observations per subject.
#' @param fixed vector of fixed effects in the following order:
#'   Intercept, main effect of A, main effect of B, and AB
#'   interaction.
#' @param amx Matrix of autocorrelation functions.
#' @param amx_wt Weights for sampling rows from \code{amx}.
#' @importFrom mgcv s
#' @export
gamm_2x2_v2 <- function(n_subj, n_obs, fixed, amx, amx_wt = NULL) {
  dat <- gen_2x2(n_subj, n_obs, fixed, TRUE, amx, amx_wt) %>%
    dplyr::mutate(subj_id = factor(subj_id),
                  tnum_r = (ts_r - (n_obs + 1) / 2) / (n_obs - 1),
                  tnum_b = (ts_b - (n_obs + 1) / 2) / (n_obs - 1))

  ## fit the models using autocorr::fit_jl
  mod_no <- mgcv::bam(Y_acn ~ AA2 * BB2 +
                        s(tnum_r, by = subj_id) +
                        s(subj_id, bs = "re") +
                        s(subj_id, AA2, bs = "re"),
                      data = dat)

  mod_rand <- mgcv::bam(Y_acr ~ AA2 * BB2 +
                          s(tnum_r, by = subj_id) +
                          s(subj_id, bs = "re") +
                          s(subj_id, AA2, bs = "re"),
                        data = dat)

  mod_block <- mgcv::bam(Y_acb ~ AA2 * BB2 +
                           s(tnum_b, by = subj_id) +
                           s(subj_id, bs = "re") +
                           s(subj_id, AA2, bs = "re"),
                         data = dat)
  
  array(c(coef(mod_no)[1:4],
          sqrt(diag(vcov(mod_no)[1:4, 1:4])),
          coef(mod_rand)[1:4],
          sqrt(diag(vcov(mod_rand)[1:4, 1:4])),
          coef(mod_block)[1:4],
          sqrt(diag(vcov(mod_block)[1:4, 1:4]))),
        dim = c(4, 2, 3),
        dimnames = list(coef = c("(Intercept)", "A", "B", "A:B"),
                        param = c("est", "stderr"),
                        model = c("mod_no", "mod_rand", "mod_block")))
}

#' Generate simulated nested categorical time-series data
#'
#' @param n_subj Number of subjects.
#' @param n_obs Number of observations per subject (multiple of 2).
#' @param fixed Fixed effects parameters (1st element = intercept, 2nd element = main effect).
#' @param downsample How much to downsample from original 1000 Hz data. (1 = 1000Hz, 2 = 500Hz, 4 = 250Hz, etc.)
#' @return Data frame with time series data.
#' @importFrom magrittr "%>%"
#' @export
gen_logit_data <- function(n_subj, n_obs, fixed, downsample = 1L) {
  dargs <- list(ivs = c(A = 2), n_item = n_obs)

  parms <- funfact::gen_pop(dargs, n_subj, var_range = c(0, .6))
  parms$item_rfx[, ] <- 0
  parms$fixed[] <- 0
  parms$err_var <- .1

  dat <- funfact::sim_norm(dargs, n_subj, parms) %>%
    tibble::as_tibble()

  dat2 <- dat %>%
    dplyr::select(-err) %>%
    dplyr::mutate(subj_id = factor(subj_id),
                  Ad = dplyr::if_else(A == "A1", -.5, .5),
                  logit = 1 / (1 + exp(-Y)),
                  ts = purrr::map(logit, oversampled_ts,
                                  downsample = downsample)) %>%
    dplyr::select(subj_id, A, Ad, Y = ts) %>%
    tidyr::unnest()
}

#' @importFrom magrittr "%>%"
#' @export
gen_2x2 <- function(n_subj, n_obs, fixed, is_acf, amx, amx_wt = NULL,
                    replace) {
  my_design <- list(ivs = c(A = 2, B = 2),
		    n_item = n_obs * 2L,
		    between_item = c("A", "B"),
		    between_subj = c("B"))
  
  parms <- funfact::gen_pop(my_design, n_subj)
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- 6

  sim_2x2(n_subj, n_obs, parms, is_acf, amx, amx_wt, replace) %>%
    funfact::with_dev_pred(c("A", "B"))
}


