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
  mx <- as.matrix(blst_studies[, c(-1, -2)])
  apply(mx, 2, quantile, probs = c(q1, q2))
}
