#' Simulate and Fit GAMM
#'
#' Simulate 2x2 data and fit GAMM/non-GAMM models
#'
#' @param n_subj Number of subjects.
#' @param n_obs Number of observations per subject.
#' @param fixed vector of fixed effects in the following order:
#'   Intercept, main effect of A, main effect of B, and AB
#'   interaction.
#' @param is_acf Is argument \code{amr} a list of residuals (TRUE) or a matrix
#'   of autocorrelation functions (FALSE)?
#' @param amr List of residuals or matrix of autocorrelation functions, with each element or row representing a different subject.
#' @param amr_wt Weights for sampling rows from \code{amx}.
#' @importFrom mgcv s
#' @export
gamm_2x2 <- function(n_subj, n_obs, fixed, is_acf, amr, amr_wt = NULL) {
  dat <- gen_2x2(n_subj, n_obs, fixed, is_acf, amr, amr_wt) %>%
    dplyr::mutate(subj_id = factor(subj_id),
                  tnum_r = (ts_r - (n_obs + 1) / 2) / (n_obs - 1),
                  tnum_b = (ts_b - (n_obs + 1) / 2) / (n_obs - 1))

  ## fit the GAMM models using mgcv::bam
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
