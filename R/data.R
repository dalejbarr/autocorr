#' Residuals from the KKL dataset (table format)
#'
#' @format A data frame with 68800 rows and 3 variables:
#' \describe{
#'   \item{subj}{subject identifier}
#'   \item{trial}{trial number (1-800)}
#'   \item{resid}{the residual}
#' }
#' @details Residuals from the KKL dataset in data frame format, after fitting a linear-mixed effects model with random intercepts and random slopes. They represent the naturally occurring autocorrelation patterns for the 86 participants in the KKL dataset, which is available through the \code{RePsychLing} package. The script for fitting the model and extracting the residuals is available at \url{https://github.com/dalejbarr/autocorr/data-raw/extract_kkl_resids.R}.
#' @source The raw KKL data can be found at \url{https://github.com/dmbates/RePsychLing}
"kkl_df"

#' Residuals from the KKL dataset (matrix format)
#'
#' @format A matrix with 86 rows and 800 columns. Each row is an individual subject, each column a trial in sequential order.
#' @details Residuals from the KKL dataset in matrix format, after fitting a linear-mixed effects model with random intercepts and random slopes. They represent the naturally occurring autocorrelation patterns for the 86 participants in the KKL dataset, which is available through the \code{RePsychLing} package. The script for fitting the model and extracting the residuals is available at \url{https://github.com/dalejbarr/autocorr/data-raw/extract_kkl_resids.R}.
#' @source The raw KKL data can be found at \url{https://github.com/dmbates/RePsychLing}
"kkl_mx"

#' Sample of study random effect variances
#'
#' @details A small sample of estimated random effect variances from
#'   13 studies, first presented in the supplementary materials to the
#'   paper by
#'   \href{http://talklab.psy.gla.ac.uk/simgen/realdata.html}{Barr,
#'   Levy, Scheepers, and Tily (2013).}
#'
#' @format A tibble with 13 rows and 6 columns:
#' \describe{
#'   \item{ID}{Numeric ID of the study (see \url{http://talklab.psy.gla.ac.uk/simgen/realdata.html}).}
#'   \item{Residual}{Estimate of the residual variance.}
#'   \item{subj_int}{By-subject random intercept variance as a proportion of the residual variance.}
#'   \item{subj_slp}{By-subject random slope variance as a proportion of the residual variance.}
#'   \item{item_int}{By-item random intercept variance as a proporition of the residual variance.}
#'   \item{item_slp}{By-item random slope variance as a proporition of the residual variance.}
#' }
"blst_studies"

#' Compute quantiles of random effect variances
#'
#' Compute quantiles of the random effect variances for the studies
#' examined by Barr, Levy, Scheepers, and Tily (2013).
#'
#' @param q1 Lower quantile.
#'
#' @param q2 Upper quantile.
#'
#' @return A matrix with the quantiles.
#' 
#' @seealso \code{\link{blst_studies}}
#'
#' @examples
#' blst_quantiles(.1, .9) # 10% and 90% cutoffs
#'
#' blst_quantiles() # 20% and 80% cutoffs
#' 
#' @export
blst_quantiles <- function(q1 = .2, q2 = .8) {
  mx <- as.matrix(autocorr::blst_studies[, c(-1, -2)])
  apply(mx, 2, quantile, probs = c(q1, q2))
}

#' Model Parameters and Residuals from Stroop ML3
#'
#' A list containing the parameter estimates and residuals from a
#' linear-mixed effects model fit to the Stroop dataset from the Many
#' Labs 3 project \insertCite{ML3}{autocorr}.
#'
#' @format A list, with the following elements:
#'
#' \describe{
#' \item{\code{fixed}}{Estimates for fixed effects, \code{(Intercept)} and
#'   main effect of congruency (\code{cong}).}
#'
#' \item{\code{covmx}}{Variance-covariance matrix for the random effects,
#'   with the random intercept and random slope variance along the
#'   diagonal, and covariances set to zero.}
#'
#' \item{\code{sigma}}{Estimate of residual variance.}
#'
#' \item{\code{param}}{A named list of 3,337 63-element vectors, each
#'   representing the model residuals for a particular participant,
#'   ordered by trial.}
#' }
#'
#' @details The model fit to the data was
#'
#' \code{lmer(latency ~ cong + (cong || session_id), stroop_ML3)}.
#'
#' @seealso \code{\link{stroop_ML3}}
#' 
#' @references
#'   \insertAllCited{}
#' 
"stroop_mod"

#' Data from the Many Labs 3 Stroop Task
#'
#' Raw data from the Stroop component of the Many Labs 3 study \insertCite{ML3}{autocorr}.
#'
#' @details Unlike the original data, the \code{latency} variable
#'   (originally called \code{trial_latency}) has been set to
#'   \code{NA} for any trials where either (1) an error occurred, or
#'   (2) the latency was greater than 10,000 milliseconds.
#' 
#' @format A tibble with 210,231 observations on 8 variables:
#'
#' \describe{
#' 
#'   \item{\code{session_id}}{Unique session (participant) identifier.}
#'
#'   \item{\code{trial}}{Trial number (0 to 62).}
#'
#'   \item{\code{ink_color}}{Color of the word.}
#'
#'   \item{\code{word}}{Identity of the word.}
#'
#'   \item{\code{congruent}}{Whether the trial was congruent or incongruent.}
#'
#'   \item{\code{cong}}{Deviation-coded predictor for \code{congruent} used in modeling.}
#'
#'   \item{\code{response}}{Color name produced by the participant.}
#'
#'   \item{\code{latency}}{Latency of the response, in milliseconds.}
#' 
#' }
#'
#' @source \url{https://osf.io/n8xa7/}
#' 
#' @examples
#' ## calculate cell means
#' with(stroop_ML3, aggregate(latency ~ congruent, stroop_ML3, mean))
#'
#' @references
#'   \insertAllCited{}
#' 
"stroop_ML3"
