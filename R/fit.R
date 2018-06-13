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
fit_jl <- function(formula, data, jsess, name = "fm1") {
  jsess$assign("form", formula)
  jsess$assign("dat", data)

  jsess$eval(paste0(name, " = fit(LinearMixedModel, form, dat)"))
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
#' @export
run_2x2 <- function(n_subj, n_obs, fixed, jsess, amx, amx_wt = NULL) {
  ## n_subj <- 48L; n_obs <- 48L; fixed <- rep(0, 4); jsess <- j
  ## rm(n_subj, n_obs, fixed, jsess)
  ## two-by-two design, everything between
  my_design <- list(ivs = c(A = 2, B = 2),
		    n_item = n_obs * 2L,
		    between_item = c("A", "B"),
		    between_subj = c("B"))

  parms <- funfact::gen_pop(my_design, n_subj)
  parms$fixed[] <- fixed
  parms$item_rfx[,] <- 0
  parms$err_var <- 6

  dat <- sim_2x2(n_subj, n_obs, parms, amx) %>%
    funfact::with_dev_pred(c("A", "B"))

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
