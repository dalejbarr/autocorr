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
