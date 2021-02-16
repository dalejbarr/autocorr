#' Simulate Practice-Effect Data and Fit All Models
#'
#' @inheritParams sim_practice
#'
#' @return A 44-element numeric vector, containing statistical results
#'   for all three models. Paramters are prefixed by \code{G} for GAMM
#'   results, by \code{L} for LMEM results, and by \code{N} for NLME
#'   results. Following the prefix, is the type of statistic, with
#'   \code{est} for parameter estimates, \code{se} for standard
#'   errors, \code{t} for t-values, and \code{p} for
#'   p-values. Finally, the suffix identifies the parameter in
#'   question, with \code{int} for the intercept, \code{wij} for the
#'   within-subject effect, and \code{bi} for the between-subject
#'   effect. There are also specific parameters for the nonlinear
#'   model, including \code{asym} for the asymptote, \code{sweep} for
#'   the sweep, and \code{lrc} for the learning rate. In cases where
#'   the non-linear model estimation fails, the error is trapped and
#'   the corresponding \code{N.xxx.xxx} values will be set to
#'   \code{NA}.
#'
#' @seealso \code{\link{sim_practice}} for information about data
#'   generation, and \code{\link{fit_practice_models}} for details
#'   about model fitting.
#'
#' @export
simfit_practice_all <- function(n_subj = 48, n_trials = 48,
				within_eff = 10,
				between_eff = 24,
				asymptote = 400,
				sweep = 400,
				learning_rate = -2,
				rslope_sd = 20,
				rasym_sd = 40,
				rsweep_sd = 40,
				rlrate_sd = 0,
				err_sd = 20,
				verbose = FALSE) {

  dat <- sim_practice(n_subj, n_trials,
		      between_eff = between_eff,
		      within_eff = within_eff,
		      asymptote = asymptote,
		      sweep = sweep,
		      learning_rate = learning_rate,
		      rslope_sd = rslope_sd, rasym_sd = rasym_sd,
		      rsweep_sd = rsweep_sd, rlrate_sd = rlrate_sd,
		      err_sd = err_sd, verbose = verbose)

  mod_gamm <- fit_practice_gamm(dat)
  mod_lmem <- fit_practice_lmem(dat)
  mod_nlme <- NULL
  tryCatch(suppressWarnings(mod_nlme <- fit_practice_nlme(dat)),
	   error = function(e) {return(NULL)})

  effs <- c("int", "wij", "bi")
  namevec <- c(paste("est", effs, sep = "."),
		     paste("se", effs, sep = "."),
		     paste("t", effs, sep = "."),
		     paste("p", effs, sep = "."))

  m_gamm_s <- summary(mod_gamm)
  gstats <- as.vector(m_gamm_s[["p.table"]])
  names(gstats) <- namevec

  m_lmem_s <- summary(mod_lmem)
  lstats <- as.vector(m_lmem_s[["p.table"]])
  names(lstats) <- namevec

  neffs <- c("wij", "asym", "bi", "sweep", "lrc")
  nstats <- rep(NA_real_, 20)
  names(nstats) <- c(paste("est", neffs, sep = "."),
		     paste("se", neffs, sep = "."),
		     paste("t", neffs, sep = "."),
		     paste("p", neffs, sep = "."))  
  if (!is.null(mod_nlme)) {
    m_nlme_s <- summary(mod_nlme)
    nstats[] <- as.vector(
      m_nlme_s[["tTable"]][,
			   c("Value", "Std.Error", "t-value", "p-value")])
  }

  c(G = gstats,
    L = lstats,
    N = nstats)
}

#' Simulate data with a practice effect
#'
#' Simulates psychological response time data containing a practice
#' effect, such that the response time decays exponentially at a rate
#' determined by \code{learning_rate} until reaching
#' \code{asymptote}. Units are milliseconds.
#'
#' @param n_subj Number of subjects.
#' @param n_trials Number of trials per subject.
#' @param within_eff Mean effect of within-subject factor, coded by
#'   \code{xi}.
#' @param between_eff Mean effect of between-subject factor, coded by
#'   \code{xj}.
#' @param asymptote Asymptotic value.
#' @param sweep The 'sweep' of the decline; i.e., the distance from
#'   the asymptote to the starting value.
#' @param learning_rate Speed of the decay.
#' @param rslope_sd By-subject standard deviation for the random slope.
#' @param rasym_sd By-subject standard deviation for the random asymptote.
#' @param rsweep_sd By-subject standard deviation for the random sweep.
#' @param rlrate_sd By-subject standard deviation for the learning rate.
#' @param err_sd Error standard deviation.
#' @param verbose Whether to return random effects in the data frame.
#'
#' @details The data contains main effects of the within-subject and
#'   between-subject factors, but no interaction. The model also
#'   contains four by-subject random effects (within-subject slope,
#'   asymptote, sweep, and learning rate). These random effects are
#'   all independent from one another. The data-generating model for
#'   response \code{y} for subject \code{i} and observation {j} is:
#'
#' \code{yij ~ (g00 + S0i) * wij + g10 + S1i + g11 * bi + 
#'       (g20 + S2i) * exp(-exp(g30 + S3i) * tij) + eij}
#'
#' where the predictors are:
#' 
#' \describe{
#'   \item{\code{g00}}{mean within-subject effect;}
#' 
#'   \item{\code{S0i}}{by-subject random slope;}
#' 
#'   \item{\code{wij}}{deviation-coded predictor for the
#'   within-subject effect;}
#' 
#'   \item{\code{g10}}{mean asymptote;}
#'
#'   \item{\code{S1i}}{random asymptote for subject \code{i};}
#' 
#'   \item{\code{g11}}{mean between-subject effect;}
#' 
#'   \item{\code{bi}}{deviation-coded predictor for the
#'   between-subject effect;}
#' 
#'   \item{\code{g20}}{mean sweep (distance from starting value
#'   to asymptote);}
#' 
#'   \item{\code{S2i}}{random sweep for subject \code{i};}
#' 
#'   \item{\code{g30}}{mean learning rate;}
#'
#'   \item{\code{S3i}}{random learning rate for subject \code{i};}
#' 
#'   \item{\code{tij}}{trial number for subject \code{i}, observation \code{j};}
#' 
#'   \item{\code{eij}}{random error for subject \code{i}, observation \code{j}.}
#' }
#' 
#' @return A data frame with some or all of the following fields,
#'   depending on the value of \code{verbose}:
#'
#' \describe{
#'   \item{\code{subj_id}}{Unique subject identifier.}
#'   \item{\code{tij}}{Trial number.}
#'   \item{\code{wij}}{Deviation-coded within-subject predictor.}
#'   \item{\code{bi}}{Deviation-coded between-subject predictor.}
#'   \item{\code{y}}{The response value.}
#'   \item{\code{S0i}}{By-subject random slope.}
#'   \item{\code{S1i}}{By-subject random asymptote.}
#'   \item{\code{S2i}}{By-subject random sweep.}
#'   \item{\code{S3i}}{By-subject random learning rate.}
#' 
#'   \item{\code{trueval}}{Response value excluding the within-subject
#'         effect/random slope.}
#' }
#'
#' @export
sim_practice <- function(n_subj, n_trials,
                      within_eff = 10, between_eff = 24,
                      asymptote = 400, sweep = 400, 
                      learning_rate = -2,
                      rslope_sd = 20,
                      rasym_sd = 40,
                      rsweep_sd = 40,
                      rlrate_sd = 0,
                      err_sd = 20,
                      verbose = FALSE) {

  genfn_ml <- function(tij, wij, bi, S0i, S1i, S2i, S3i,
                       g00, g10, g11, g20, g30) {
    (g00 + S0i) * wij + g10 + S1i + g11 * bi + 
      (g20 + S2i) * exp(-exp(g30 + S3i) * tij)
  }

  trial <- seq_len(n_trials)

  dat <- data.frame(
    subj_id = factor(rep(sprintf("S%02d", seq_len(n_subj)),
                         each = length(trial))),
    tij = rep(trial, times = n_subj),
    wij = unlist(replicate(n_subj, {
      sample(rep(c(-.5,.5), each = length(trial) / 2))
    }, simplify=FALSE)), # randomize trial order
    bi = rep(sample(rep(c(-.5, .5), each = n_subj / 2)),
             each = length(trial)),
    S0i = rep(rnorm(n_subj, sd = rslope_sd), each = length(trial)),
    S1i = rep(rnorm(n_subj, sd = rasym_sd), each = length(trial)),
    S2i = rep(rnorm(n_subj, sd = rsweep_sd), each = length(trial)),
    S3i = rep(rnorm(n_subj, sd = rlrate_sd), each = length(trial)))

  dat[["trueval"]] <- with(
    dat,
    mapply(genfn_ml,
           tij, 0, bi, 0, S1i, S2i, S3i,
           g00 = within_eff, g11 = between_eff,
           g10 = asymptote, g20 = sweep, g30 = learning_rate))

  dat[["y"]] <- with(
    dat,
    mapply(genfn_ml,
           tij, wij, bi, S0i, S1i, S2i, S3i,
           g00 = within_eff, g11 = between_eff,
           g10 = asymptote, g20 = sweep, g30 = learning_rate)) +
    rnorm(n_subj * length(trial), sd = 20)

  if (verbose) {
    dat
  } else {
    dat[, c("subj_id", "tij", "wij", "bi", "y")]
  }
}


#' Fit Models to Simulated Practice Data
#'
#' @name fit_practice_models
#'
#' @param dat Data frame, the result of \code{\link{sim_practice}}.
#' 
#' @param start Vector of five starting values for
#'   \code{\link[nlme]{nlme}}, for the values (1) within-effect; (2)
#'   asymptote; (3) between-effect; (4) sweep; and (5) learning rate.
#' 
#' @details \code{fit_practice_nlme} fits a Non-Linear Mixed-Effects
#'   (NLME) model using \code{\link[nlme]{nlme}};
#'   \code{fit_practice_gamm} fits a Generalized Additive
#'   Mixed-Effects Model (GAMM) using \code{\link[mgcv]{bam}} with a
#'   by-subject factor smooth for trial; \code{fit_practice_lmem} fits
#'   a Linear Mixed-Effects Model (LMEM) also using
#'   \code{\link[mgcv]{bam}} for comparability to the GAMM, but as it
#'   includes no wiggly terms, it is functionally equivalent to a
#'   model fit using the \code{lme4} package. For simplicity, all
#'   models assume a constant learning rate across subjects.
NULL

#' @rdname fit_practice_models
#' 
#' @export
fit_practice_nlme <- function(dat,
                    start_vals = c(10, 400, 24, 400, -2)) {
  
  nlme::nlme(y ~ b0 * wij + Asym + sweep * exp(-exp(lrc) * tij),
                      groups = ~ subj_id,
                      dat, fixed = list(b0 ~ 1, Asym ~ bi, sweep ~ 1, lrc ~ 1),
                      random = nlme::pdDiag(b0 + Asym + sweep ~ 1),
                      start = start_vals)
}

#' @rdname fit_practice_models
#'
#' @param dat Data frame, the result of \code{\link{sim_practice}}.
#' 
#' @details 
#' 
#' @export
fit_practice_gamm <- function(dat) {
  mgcv::bam(y ~ wij + bi +
              s(tij, bs = "tp") + # common smooth
              s(wij, subj_id, bs = "re") + # random slope
              s(tij, subj_id, bs = "fs"), # factor smooth (time-varying icept)
            data = dat)
}

#' @rdname fit_practice_models
#'
#' @export
fit_practice_lmem <- function(dat) {
  mgcv::bam(y ~ wij + bi +
              s(subj_id, bs = "re") + # random intercept
              s(wij, subj_id, bs = "re"), # random slope
            data = dat)
}

