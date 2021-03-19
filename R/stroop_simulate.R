#' Simulate Stroop Factorial Data
#'
#' Simulate data from a mixed factorial design based on the Many Labs
#' 3 (ML3) \insertCite{ML3}{autocorr} Stroop data.
#'
#' @param n_subj Number of subjects.
#'
#' @param B1 Fixed effect of congruency.
#'
#' @param B2 Fixed effect of the between-subjects factor.
#'
#' @param rivar Random intercept variance.
#'
#' @param rsvar Random slope variance.
#'
#' @param verbose Whether to return error components (random effects and residuals).
#' 
#' @details The ML3 data was from a simple design including only the
#'   within-factor of congruency. The simulated data includes an
#'   additional between-subjects factor. The simulation parameters are
#'   taken from \code{\link{stroop_mod}}, with a random subset of
#'   residuals grafted onto the fitted values for each subject. 
#'
#' The simulated response latency for observation \code{j} of subject
#' \code{i} is given by the following formula:
#'
#' \code{latency_ij = B0 + b_0i + B2 * B_i + (B1 + b_1i) * W_ij + e_ij}
#'
#' where \code{B0} is \code{stroop_mod$fixed["(Intercept)"]},
#' \code{b_0i} and \code{b_1i} are the random intercept and random
#' slope for subject \code{i}, \code{W_ij} and \code{B_i} are
#' deviation-coded predictors for the within- and between- subject
#' factors respectively, \code{B1} and \code{B2} are the raw within- and
#' between- subject effect sizes, and \code{e_ij} is the \code{j}th
#' residual for subject \code{i}.
#'
#' Note that whereas the original dataset had 63 observations per
#' subject (with 2/3 in the incongruent condition and 1/3 in the
#' congruent condition), the current dataset had 1/2 congruent and 1/2
#' incongruent. Also, only the first 20 of the 63 residuals were used,
#' to eliminate discontinuities introduced at the start of each
#' 21-trial block.
#'
#' @return A data frame with \code{n_subj * 20} simulated observations
#'   on 7 (or 10) variables, depending on \code{verbose}:
#'
#' \describe{
#'
#'   \item{\code{session_id}}{A factor with \code{session_id} values
#'   corresponding to the subjects whose residuals were sampled from
#'   \code{\link{stroop_mod[["resid"]]}}}.
#'
#'   \item{\code{trial}}{An integer specifying the trial number (1 to 20).}
#'
#'   \item{\code{congruency}}{Level of within-subject factor (\code{congruent} or \code{incongruent}).}
#'
#'   \item{\code{group}}{Level of between-subject factor (\code{G1} or \code{G2}).}
#'
#'   \item{\code{W_ij}}{Deviation-coded predictor for factor \code{congruency}.}
#'
#'   \item{\code{B_i}}{Deviation-coded predictor for factor \code{group}.}
#'
#'   \item{\code{Y_ij}}{Simulated response latency in milliseconds.}
#'
#'   \item{\code{b_0i}}{Random effect for subject i (verbose mode only).}
#'
#'   \item{\code{b_1i}}{Random slope for subject i (verbose mode only).}
#'
#'   \item{\code{e_ij}}{Residual for subject i, observation j (verbose mode only).}
#' }
#'
#' @examples
#' ## use parameters from Many Labs 3 model fit
#' dat <- simulate_stroop(48, A = stroop_mod$fixed["cong"])
#'
#' mod <- lme4::lmer(Y_ij ~ W_ij + B_i + (W_ij || session_id), dat)
#' summary(mod)
#' 
#' @export
simulate_stroop <- function(n_subj, B1 = 0, B2 = 0,
                            rivar = diag(autocorr::stroop_mod$covmx)[1],
                            rsvar = diag(autocorr::stroop_mod$covmx)[2],
                            verbose = FALSE) {
  sids <- sample(names(autocorr::stroop_mod[["resid"]]), n_subj)
  sids <- sids[order(as.integer(sids))]
  
  srfx <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = autocorr::stroop_mod[["covmx"]])

  dat <- data.frame(session_id = rep(factor(sids), each = 20L),
                    trial = rep(1:20, times = n_subj),
                    b_0i = rep(srfx[, 1], each = 20L),
                    b_1i = rep(srfx[, 2], each = 20L),
                    congruency = factor(
                      unlist(replicate(n_subj,
                                       sample(rep(c("incongruent", "congruent"),
                                                  each = 10L)),
                                       simplify = FALSE))),
                    group = factor(rep(sample(rep(c("G1", "G2"), each = n_subj / 2L)),
                                   each = 20L)))
  
  dat[["W_ij"]] <- ifelse(dat[["congruency"]] == "incongruent", -.5, .5)
  dat[["B_i"]] <- ifelse(dat[["group"]] == "G1", -.5, .5)
  dat[["e_ij"]] <- unlist(lapply(sids, function(nx) {
    autocorr::stroop_mod[["resid"]][[nx]][1:20]
  }))
  dat[["Y_ij"]] <- autocorr::stroop_mod[["fixed"]]["(Intercept)"] + dat[["b_0i"]] +
    B2 * dat[["B_i"]] +
    (B1 + dat[["b_1i"]]) * dat[["W_ij"]] +
    dat[["e_ij"]]

  keep_cols <- c("session_id", "trial", "congruency", "group", "W_ij", "B_i", "Y_ij")

  if (verbose) {
    dat[, c(keep_cols, c("b_0i", "b_1i", "e_ij"))]
  } else {
    dat[, keep_cols]
  }
}

#' Fit GAMM and LMEM Models to Simulated Stroop Data
#'
#' Uses \code{\link[mgcv]{bam}} to fit two Generalized Additive
#' Mixed-Effects Models (GAMM) and a Linear Mixed-Effects Model (LMEM)
#' to simulated Stroop data.
#'
#' @details The models fit to the data are:
#'
#' @section GAMM:
#'
#' with penalized factor smooths:
#' 
#' \code{bam(Y_ij ~ W_ij + B_i +
#'       s(trial, bs = "tp") +                    # common smooth
#'       s(session_id, trial, bs = "fs", m = 1) + # factor smooth
#'       s(W_ij, session_id, bs = "re"),           # random slope
#'       data = dat)}
#'
#' and unpenalized factor smooths:
#'
#' \code{bam(Y_ij ~ W_ij + B_i +
#'       s(trial, bs = "tp") +            # common smooth
#'       s(session_id, trial, bs = "fs")  # factor smooth
#'       s(W_ij, session_id, bs = "re"),  # random slope
#'       data = dat)}
#' 
#' @section LMEM:
#' 
#' \code{bam(Y_ij ~ W_ij + B_i +
#'       s(session_id, bs = "re") +      # random intercept
#'       s(W_ij, session_id, bs = "re"), # random slope
#'       data = dat)}
#' 
#' @param dat A data.frame with simulated stroop data, the result of
#'   \code{\link{simulate_stroop}}.
#'
#' @return A named 27-element vector with the following statistical results:
#' 
#' \describe{
#'
#'   \item{\code{Gp.e.B0}}{Estimated fixed intercept of the penalized GAMM model.}
#'
#'   \item{\code{Gp.e.B1}}{Estimated fixed effect of within-subject
#'   factor of the penalized GAMM model.}
#'
#'   \item{\code{Gp.e.B2}}{Estimated fixed effect of between-subject
#'   factor of the penalized GAMM model.}
#'
#'   \item{\code{Gp.se.B0}}{Standard error for the intercept of the penalized GAMM model.}
#'
#'   \item{\code{Gp.se.B1}}{Standard error for within-effect of the penalized GAMM model.}
#'
#'   \item{\code{Gp.se.B2}}{Standard error for between-effect of the penalized GAMM model.}
#'
#'   \item{\code{Gp.p.B0}}{P-value for the intercept of the penalized GAMM model.}
#'
#'   \item{\code{Gp.p.B1}}{P-value for within-effect of the penalized GAMM model.}
#'
#'   \item{\code{Gp.p.B2}}{P-value for between-effect of the penalized GAMM model.}
#'
#'   \item{\code{Gu.e.B0}}{Estimated fixed intercept of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.e.B1}}{Estimated fixed effect of within-subject
#'   factor of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.e.B2}}{Estimated fixed effect of between-subject
#'   factor of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.se.B0}}{Standard error for the intercept of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.se.B1}}{Standard error for within-effect of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.se.B2}}{Standard error for between-effect of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.p.B0}}{P-value for the intercept of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.p.B1}}{P-value for within-effect of the unpenalized GAMM model.}
#'
#'   \item{\code{Gu.p.B2}}{P-value for between-effect of the unpenalized GAMM model.}
#'
#'   \item{\code{Lm.e.B0}}{Estimated fixed intercept of the LMEM model.}
#'
#'   \item{\code{Lm.e.B1}}{Estimated fixed effect of within-subject
#'   factor of the LMEM model.}
#'
#'   \item{\code{Lm.e.B2}}{Estimated fixed effect of between-subject
#'   factor of the LMEM model.}
#'
#'   \item{\code{Lm.se.B0}}{Standard error for the intercept of the LMEM model.}
#'
#'   \item{\code{Lm.se.B1}}{Standard error for within-effect of the LMEM model.}
#'
#'   \item{\code{Lm.se.B2}}{Standard error for between-effect of the LMEM model.}
#'
#'   \item{\code{Lm.p.B0}}{P-value for the intercept of the LMEM model.}
#'
#'   \item{\code{Lm.p.B1}}{P-value for within-effect of the LMEM model.}
#'
#'   \item{\code{Lm.p.B2}}{P-value for between-effect of the LMEM model.}
#'
#' }
#'
#' @examples
#' \donttest{
#'   fit_stroop(simulate_stroop(20))
#' }
#' 
#' @export
fit_stroop <- function(dat) {
  
  mod_gam_p <- mgcv::bam(Y_ij ~ W_ij + B_i +
                           s(trial, bs = "tp") +
                           s(session_id, trial, bs = "fs", m = 1) +
                         s(W_ij, session_id, bs = "re"), data = dat)
  
  mod_gam_u <- mgcv::bam(Y_ij ~ W_ij + B_i +
                           s(trial, bs = "tp") +
                           s(session_id, trial, bs = "fs") +
                           s(W_ij, session_id, bs = "re"), data = dat)

  mod_lmm <- mgcv::bam(Y_ij ~ W_ij + B_i +
                         s(session_id, bs = "re") +
                         s(W_ij, session_id, bs = "re"), data = dat)

  mod_lmm_s <- summary(mod_lmm)
  mod_gam_p_s <- summary(mod_gam_p)
  mod_gam_u_s <- summary(mod_gam_u)
  
  vv_lmm <- c(coef(mod_lmm)[1:3],
              sqrt(diag(vcov(mod_lmm)[1:3, 1:3])),
              mod_lmm_s[["p.table"]][, "Pr(>|t|)"])

  names(vv_lmm) <- c("e.B0", "e.B1", "e.B2",
                     "se.B0", "se.B1", "se.B2",
                     "p.B0", "p.B1", "p.B2")

  vv_gam_p <- c(coef(mod_gam_p)[1:3],
    sqrt(diag(vcov(mod_gam_p)[1:3, 1:3])),
    mod_gam_p_s[["p.table"]][, "Pr(>|t|)"])

  vv_gam_u <- c(coef(mod_gam_u)[1:3],
    sqrt(diag(vcov(mod_gam_u)[1:3, 1:3])),
    mod_gam_u_s[["p.table"]][, "Pr(>|t|)"])
  
  names(vv_gam_p) <- names(vv_lmm)
  names(vv_gam_u) <- names(vv_lmm)

  c(Gp = vv_gam_p,
    Gu = vv_gam_u,
    Lm = vv_lmm)
}
