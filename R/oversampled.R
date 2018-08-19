#' Generate single oversampled categorical time series
#'
#' @param p_targ Probability of gazing at the target.
#' @param prob Probability of a decision at each frame; 1 / prob gives mean lag between decisions.
#' @param downsample 1L for 1000 Hz; 2L for 500 Hz; 4L for 250 Hz, etc.
#' @param sacc_overhead Saccadic overhead.
#' @details Decision lags are sampled from a geometric distribution with probability \code{prob}, with each decision offset from the last by \code{sacc_overhead}.
#' @return Oversampled time series, represented by a vector of zeroes and ones 
#' @export
oversampled_ts <- function(p_targ,
			   prob = .005,
			   downsample = 1L,
			   sacc_overhead = 200L) {
  if (!is.integer(downsample))
    stop("'downsample' must be an integer (e.g., 1L, 2L, 3L)")
  n_frames_per_trial = 1000L

  npts <- ceiling(n_frames_per_trial / sacc_overhead) + 1L

  vec <- rgeom(npts, prob) + c(0L, rep(sacc_overhead, npts - 1L))

  n_decisions <- min(which(cumsum(vec) > 1000L))

  ## truncate vec and equate to n_frames_per_trial
  vec2 <- vec[seq_len(n_decisions)]
  vec2[length(vec2)] <- vec2[length(vec2)] - (sum(vec2) - n_frames_per_trial)

  ## decisions
  pog <- sample(c(1L, 0L), n_decisions, TRUE, c(p_targ, 1 - p_targ))

  frames <- rep(pog, vec2)
  frames[seq(1L, n_frames_per_trial, downsample)]
}
