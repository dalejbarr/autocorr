#' Simulate Stroop Factorial Data
#'
#' Simulate data from a mixed factorial design based on the Many Labs
#' 3 (ML3) \insertCite{ML3}{autocorr} Stroop data.
#'
#' @param n_subj Number of subjects.
#'
#' @param A Fixed effect of congruency.
#'
#' @param B Fixed effect of the between-subjects factor.
#'
#' @param rivar Random intercept variance.
#'
#' @param rsvar Random slope variance.
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
#' \code{latency_ij = int + sri_i + B * B_c + (A + srs_i) * A_c + e_ij}
#'
#' where \code{int} is \code{stroop_mod$fixed["(Intercept)"]},
#' \code{sri_i} and \code{srs_i} are the random intercept and random
#' slope for subject \code{i}, \code{A_c} and \code{B_c} are
#' deviation-coded predictors for the within- and between- subject
#' factors respectively, \code{A} and \code{B} are the raw within- and
#' between- subject effect sizes, and \code{e_ij} is the \code{j}th
#' residual for subject \code{i}.
#'
#' Note that whereas the original dataset had 63 observations per
#' subject (with 2/3 in the incongruent condition and 1/3 in the
#' congruent condition), the current dataset had 1/2 congruent and 1/2
#' incongruent. To make this work, the 63rd residual was lopped off
#' for each participant.
#'
#' @return
#' A data frame with \code{n_subj * 62} simulated observations on
#' \describe{
#'
#'   \item{\code{session_id}}{A factor with \code{session_id} values
#'   corresponding to the subjects whose residuals were sampled from
#'   \code{\link{stroop_mod[["resid"]]}}}.
#'
#'   \item{\code{trial}}{An integer specifying the trial number (0 to 62).}
#'
#'   \item{\code{A}}{Level of within-subject factor (\code{A1} or \code{A2}).}
#'
#'   \item{\code{B}}{Level of between-subject factor (\code{B1} or \code{B2}).}
#'
#'   \item{\code{A_c}}{Deviation-coded predictor for factor \code{A}.}
#'
#'   \item{\code{B_c}}{Deviation-coded predictor for factor \code{B}.}
#'
#'   \item{\code{latency}}{Simulated response latency in milliseconds.}
#' }
#'
#' @examples
#' ## use parameters from Many Labs 3 model fit
#' dat <- simulate_stroop(48, A = stroop_mod$fixed["cong"])
#'
#' mod <- lme4::lmer(latency ~ A_c + B_c + (A_c || session_id), dat)
#' summary(mod)
#' 
#' @export
simulate_stroop <- function(n_subj, A = 0, B = 0,
                            rivar = diag(autocorr::stroop_mod$covmx)[1],
                            rsvar = diag(autocorr::stroop_mod$covmx)[2]) {
  sids <- sample(names(autocorr::stroop_mod[["resid"]]), n_subj)
  sids <- sids[order(as.integer(sids))]
  
  srfx <- MASS::mvrnorm(n_subj, mu = c(0, 0), Sigma = autocorr::stroop_mod[["covmx"]])

  dat <- data.frame(session_id = rep(factor(sids), each = 62L),
                    trial = rep(0:61, times = n_subj),
                    sri = rep(srfx[, 1], each = 62L),
                    srs = rep(srfx[, 2], each = 62L),
                    A = factor(unlist(replicate(n_subj, sample(rep(c("A1", "A2"), each = 31L)),
                                         simplify = FALSE))),
                    B = factor(rep(sample(rep(c("B1", "B2"), each = n_subj / 2L)), each = 62L)))
  
  dat[["A_c"]] <- ifelse(dat[["A"]] == "A1", -.5, .5)
  dat[["B_c"]] <- ifelse(dat[["B"]] == "B1", -.5, .5)
  dat[["err"]] <- unlist(lapply(sids, function(nx) {
    autocorr::stroop_mod[["resid"]][[nx]][1:62]
  }))
  dat[["latency"]] <- autocorr::stroop_mod[["fixed"]]["(Intercept)"] + dat[["sri"]] +
    B * dat[["B_c"]] +
    (A + dat[["srs"]]) * dat[["A_c"]] +
    dat[["err"]]
  dat[, c("session_id", "trial", "A", "B", "A_c", "B_c", "latency")]
}

#' Fit GAMM and LMEM Models to Simulated Stroop Data
#'
#' Uses \code{\link[mgcv]{bam}} to fit a Generalized Additive
#' Mixed-Effects Model (GAMM) and a Linear Mixed-Effects Model (LMEM)
#' to simulated Stroop data.
#'
#' @details The models fit to the data are:
#'
#' @section GAMM:
#' \code{bam(latency ~ A_c + B_c +
#'       s(trial, bs = "tp") +                    # common smooth
#'       s(session_id, trial, bs = "fs", m = 1) + # factor smooth
#'       s(A_c, session_id, bs = "re"),           # random slope
#'       data = dat)}
#'
#' @section LMEM:
#' \code{bam(latency ~ A_c + B_c +
#'       s(session_id, bs = "re") +  # random intercept
#'       s(A_c, session_id, bs = "re"), # random slope
#'       data = dat)}
#' 
#' @param dat A data.frame with simulated stroop data, the result of
#'   \code{\link{simulate_stroop}}.
#'
#' @return A named vector with the following statistical results:
#' 
#' \describe{
#'
#'   \item{\code{G.e.int}}{Estimated fixed intercept of the GAMM model.}
#'
#'   \item{\code{G.e.A}}{Estimated fixed effect of A of the GAMM model.}
#'
#'   \item{\code{G.e.B}}{Estimated fixed effect of B of the GAMM model.}
#'
#'   \item{\code{G.se.int}}{Standard error for the intercept of the GAMM model.}
#'
#'   \item{\code{G.se.A}}{Standard error for effect A of the GAMM model.}
#'
#'   \item{\code{G.se.B}}{Standard error for effect B of the GAMM model.}
#'
#'   \item{\code{G.p.int}}{P-value for the intercept of the GAMM model.}
#'
#'   \item{\code{G.p.A}}{P-value for A of the GAMM model.}
#'
#'   \item{\code{G.p.B}}{P-value for B of the GAMM model.}
#'
#'   \item{\code{L.e.int}}{Estimated fixed intercept of the LMEM model.}
#'
#'   \item{\code{L.e.A}}{Estimated fixed effect of A of the LMEM model.}
#'
#'   \item{\code{L.e.B}}{Estimated fixed effect of B of the LMEM model.}
#'
#'   \item{\code{L.se.int}}{Standard error for the intercept of the LMEM model.}
#'
#'   \item{\code{L.se.A}}{Standard error for effect A of the LMEM model.}
#'
#'   \item{\code{L.se.B}}{Standard error for effect B of the LMEM model.}
#'
#'   \item{\code{L.p.int}}{P-value for the intercept of the LMEM model.}
#'
#'   \item{\code{L.p.A}}{P-value for A of the LMEM model.}
#'
#'   \item{\code{L.p.B}}{P-value for B of the LMEM model.}
#'
#' }
#' 
#' @export
fit_stroop <- function(dat) {
  mod_gam <- mgcv::bam(latency ~ A_c + B_c +
                         s(trial, bs = "tp") +
                         s(session_id, trial, bs = "fs", m = 1) +
                         s(A_c, session_id, bs = "re"), data = dat)

  mod_lmm <- mgcv::bam(latency ~ A_c + B_c +
                         s(session_id, bs = "re") +
                         s(A_c, session_id, bs = "re"), data = dat)

  mod_gam_s <- summary(mod_gam)
  mod_lmm_s <- summary(mod_lmm)

  vv_gam <- c(coef(mod_gam)[1:3],
    sqrt(diag(vcov(mod_gam)[1:3, 1:3])),
    mod_gam_s[["p.table"]][, "Pr(>|t|)"])

  names(vv_gam) <- c("e.int", "e.A", "e.B",
                     "se.int", "se.A", "se.B",
                     "p.int", "p.A", "p.B")

  vv_lmm <- c(coef(mod_lmm)[1:3],
    sqrt(diag(vcov(mod_lmm)[1:3, 1:3])),
    mod_lmm_s[["p.table"]][, "Pr(>|t|)"])

  names(vv_lmm) <- c("e.int", "e.A", "e.B",
                     "se.int", "se.A", "se.B",
                     "p.int", "p.A", "p.B")

  c(G = vv_gam,
    L = vv_lmm)
}
