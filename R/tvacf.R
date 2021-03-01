#' Time-Varying Autocorrelation Function
#'
#' @param x A time series.
#'
#' @details Slides a moving window along vector \code{x} with
#'   specified \code{width} and calulates the autocorrelation function
#'   at each time slice.
#'
#' @return A matrix, the columns vectors of which represent the
#'   results of \code{\link[stats]{acf}} applied to each time slice,
#'   including lags from 1 to \code{width} - 1.
#'
#' @examples
#' acseries <- stat_gp(1, 2)$GP[1:50]
#' tvacf(acseries - mean(acseries), width = 10)
#'  
#' @export
tvacf <- function(x, width = 7L) {
  if (width >= length(x)) {
    stop("'width' must be smaller than the length of 'x'")
  }

  last_x <- length(x) - width
  sapply(seq_len(last_x), function(.x) {
    stats::acf(x[.x:(.x + width - 1L)], plot = FALSE)[[1]][-1, , 1]
  })
}

tvacf_sig <- function(x, width = 7L, nmc = 1000L) {
  orig <- tvacf(x, width)

  pvals <- sapply(1:(length(x) - width), function(t0) {
    chunk <- x[t0:(t0 + width - 1L)]
    nhd <- replicate(nmc, {
      stats::acf(sample(chunk), plot = FALSE)[[1]][-1, , 1]
    })

    apply(cbind(orig[, t0], nhd), 1, function(.x) {
      sum(abs(.x) >= abs(.x)[1]) / length(.x)
    })
  })

  pvals * sign(orig)
}

tvacf_clusters <- function(x, width = 7L,
                           nmc_outer = 4L,
                           nmc_inner = 1000L,
                           alpha_outer = .05,
                           alpha_inner = .05,
                           center = TRUE) {
  v <- if (center) {
         x - mean(x)
       } else {
         x
       }
  
  orig_acf <- tvacf(v, width)
  orig_p  <- tvacf_sig(v, width, nmc_inner)

  nhd <- replicate(nmc_outer, {
    tvacf_sig(sample(v), width, nmc_inner)
  })

  list(orig_acf = orig_acf,
       orig_p = orig_p,
       nhd = nhd)
}

## cluster mass statistics
detect_clusters <- function(x, alpha = .05) {
  runs <- rle(as.integer((abs(x) < alpha) * sign(x)))
  clusters <- which(abs(runs$values) == 1L)
  if (length(clusters) == 0L) {
    data.frame(t0 = NA_integer_, t1 = NA_integer_,
               sign = NA_integer_, cms = 0)
  } else {
    cluster_ix <- sapply(clusters, function(.x) {
      if (.x == 1L) {
        c(1L, runs$lengths[.x])
      } else {
        c(sum(runs$lengths[seq_len(.x) - 1L]) + 1L,
          sum(runs$lengths[seq_len(.x) - 1L]) + runs$lengths[.x])
      }
    })
    vals <- sapply(seq_len(dim(cluster_ix)[2]), function(.cx) {
      .x <- cluster_ix[, .cx]
      sum( -log(abs(x[.x[1]:.x[2]])) )
    })
    data.frame(t0 = cluster_ix[1, ],
               t1 = cluster_ix[2, ],
               sign = sign(x[cluster_ix[1, ]]),
               cms = vals)
  }
}
