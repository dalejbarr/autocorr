## residual function with exponential decay
expdecay <- function(n_subj, n_obs,
                     lambda = runif(1, log(2), log(15))) {
  ## subjects are normally distributed around lambda
  offset <- rnorm(n_subj, sd = .8)

  lambda_i <- exp(lambda + offset)

  res <- lapply(lambda_i,
                function(.x) exp(-.x * seq(0, 1, length.out = n_obs)))

  ## normalize
  nres <- lapply(res,
                 function(.x) (.x - mean(.x)) / sd(.x))

  ## add small amount of noise
  nres_noise <- lapply(nres,
                       function(.x) .x + rnorm(n_obs, sd = sqrt(.1)))

  ## normalize again
  lapply(nres_noise, function(.x) (.x - mean(.x)) / sd(.x))
}

