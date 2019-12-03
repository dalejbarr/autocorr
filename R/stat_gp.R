#' Generate autocorrelation function
#'
#' @details Autocorrelation function for a Gaussian process, defined as \code{sigma^2 * exp(-t^2 / (2 * gamma^2))}. To convert covariances to correlations, use \code{cov_fn(t, sigma, gamma) / sigma^2}.
#' 
#' @param t Vector of time lags.
#' @param sigma Standard deviation for the Gaussian process.
#' @param gamma Parameter co-determining (with sigma) the length of the autocorrelation.
#' @return Vector of covariances. The duration of the autocorrelation across lags is determined by the ratio of sigma and gamma.
#' @examples
#' sig <- 4
#' plot(cov_fn(1:20, sigma = sig) / sig^2,
#'      type = 'b', ylim = c(0, 1), ylab = 'correlation', xlab = 'lag')
#' @export
cov_fn <- function(t, sigma = 4, gamma = 3) {
  sigma^2 * exp(-t^2 / (2 * gamma^2))
}


#' Generate Gaussian time-series with autocorrelation
#'
#' @param sigma Standard deviation for the gaussian.
#' @param gamma The correlation length.
#' @param dt_GP Desired time step for the series.
#' @details TODO
#' @return TODO
#' @export
stat_gp <- function(sigma = 4, gamma = 3, dt_GP = 0.8) {
  ## utility function
  cumsimpson <- function(x, y) {
    ry <- 1L
    cy <- length(y)
    
    dx <- diff(x)
    dy <- diff(y)
    dx1 <- dx[-length(dx)]
    dx2 <- dx[-1L]
    dxs <- dx1 + dx2
    dy1 <- dy[-length(dy)]
    dy2 <- dy[-1L]
    a <- (dy2 / (dx2 * dxs) - dy1 / (dx2 * dxs)) / 3
    b <- (dy2 * dx1 / (dx2 * dxs) + dy1 * dx2 / (dx1 * dxs)) / 2
    c <- y[-c(1L, length(y))]
    
    i1 <- ((a * dx1 - b) * dx1 + c) * dx1 # left half integral
    i2 <- ((a * dx2 + b) * dx2 + c) * dx2 # right half integral
    
    z <- vector("numeric", length(y))
    z[seq(2, length(z) - 1, 2)] <- i1[seq(1, length(i1), 2)]
    z[seq(3, length(z), 2)] <- i2[seq(1, length(i2), 2)]
    z[length(z)] <- i2[length(i2)]
    cumsum(z)
  }

  ## Generates a GP by choosing dt in the spectral representation. This is the
  ## latest version of the implementation that unties variables for the
  ## spectrum computation from those for generating the GP.
    
  ## compute spectrum and spectral cut-off frequency
  ## cov_fn.param <- formals(cov_fn)
  Nt_cor <- 100
  dt <- gamma / Nt_cor
  
  eps.cov <- 1e-4
  ## find the root of the function, store in Tmax$x
  Tmax <- pracma::fzero(function(t) cov_fn(t, sigma, gamma) - eps.cov,
                        5 * gamma, tol = 1e-8)
  
  Tmax <- ceiling(Tmax$x / dt) * dt
  Nt <- Tmax / dt
  ind_0 <- Nt + 1
  k <- -Nt:(Nt - 1)
  tcov <- k * dt
  cov.t <- cov_fn(tcov, sigma, gamma)
  
  freq <- k / (2 * Tmax)
  w <- 2 * pi * freq
  dw <- 2 * pi / (2 * Tmax) ## equation (43)
  S <- abs(pracma::fftshift(fft(cov.t)) * Tmax / Nt) / (2 * pi)
  Int_S <- cumsimpson(w,S)
  
  ## Test for convergence of spectral integral
  rel_error_Int_S = (Int_S[(2*Nt - 10L):(2*Nt)] -
                       Int_S[(2*Nt - 11L):(2*Nt - 1L)]) /
    Int_S[(2*Nt - 10L):(2*Nt)]
  eps_Int_S <- 1e-3
  if (sum(abs(rel_error_Int_S) > eps_Int_S) > 0) {
    stop("w-range for spectrum too small. Increase Tmax. (eps_Int_S = ",
         eps_Int_S, ").")
  }
  
  eps_cutoff <- 1e-3
  
  ind_wc <- which(Int_S > (1 - eps_cutoff) * Int_S[2*Nt])[1]
  Nw <- ind_wc - ind_0
  wc <- Nw * dw

  ## alter this process
  ## generate stochastic process with dt close to intial guess di_GP
  Nw <- 128
  dw2 <- wc/Nw
  T0 <- 2 * pi / dw2 ## total length of the GP, controlled by Mw
  Mw <- floor(T0 / dt_GP)
  dt <- T0 / Mw ## fixed dt

  ## critical test that everything makes sense
  if (dt > pi / wc)
    stop('Reduce dt to be consistent with cut-off frequency wc')
  
  ## recompute spectrum with new spectral discretisation
  dts <- 2 * pi / (2 * wc)
  cov_t <- cov_fn((-Nw:(Nw - 1)) * dts, sigma, gamma)
  S2 <- abs((fft(cov_t) * dts)) / (2 * pi)
  w2 <- (-Nw:(Nw - 1)) * dw2
  
  A <- vector("numeric", Mw)
  A[1:Nw] <- sqrt(2 * S2[1:Nw] * dw2)
  A[1] <- 0
  
  phi <- vector("numeric", Mw)
  phi[1:Nw] <- 2 * pi * runif(Nw)
  B <- sqrt(2) * A * exp(complex(real = 0, imaginary = 1) * phi)
  
  GP <- Mw * Re(fft(B, inverse = TRUE) / length(B))
  return(list("GP" = GP, "dt" = dt))
}

