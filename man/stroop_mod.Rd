% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{stroop_mod}
\alias{stroop_mod}
\title{Model Parameters and Residuals from Stroop ML3}
\format{
A list, with the following elements:

\describe{
\item{\code{fixed}}{Estimates for fixed effects, \code{(Intercept)} and
main effect of congruency (\code{cong}).}

\item{\code{covmx}}{Variance-covariance matrix for the random effects,
with the random intercept and random slope variance along the
diagonal, and covariances set to zero.}

\item{\code{sigma}}{Estimate of residual variance.}

\item{\code{param}}{A named list of 3,337 63-element vectors, each
representing the model residuals for a particular participant,
ordered by trial.}
}
}
\usage{
stroop_mod
}
\description{
A list containing the parameter estimates and residuals from a
linear-mixed effects model fit to the Stroop dataset from the Many
Labs 3 project \insertCite{ML3}{autocorr}.
}
\details{
The model fit to the data was

\code{lmer(latency ~ cong + (cong || session_id), stroop_ML3)}.
}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{stroop_ML3}}
}
\keyword{datasets}
