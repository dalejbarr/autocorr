#' Simulate 2x2 data with autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate. Must be a positive
#'   integer that is a multiple of 2.
#' 
#' @param n_items Number of items/observations per subject. Must be a positive
#'   even integer.
#'
#' @param int Intercept.
#'
#' @param A Main effect of A (within-subject factor).
#'
#' @param B Main effect of B (between-subject factor).
#'
#' @param AB AB interaction effect.
#'
#' @param item_rint Item random intercept variance.
#'
#' @param item_rslp Item random slope variance (for factor B).
#'
#' @param item_rcorr Item random correlation.
#'
#' @param subj_rint Subject random intercept variance.
#'
#' @param subj_rslp Subject random slope variance (for factor A).
#'
#' @param subj_rcorr Subject random correlation.
#'
#' @param version How to generate residuals: either an integer
#'   representing the Scenario number (see \code{\link{errsim}}) or
#'   the name of a user-defined function.
#'
#' @param verbose Whether the data frame should GLM components.
#' 
#' @param extra_args A list of extra arguments to be passed to a
#'   user-defined function to generate residuals.
#'
#' @return A data frame with \code{n_subj * n_items} rows and either 9
#'   or 12 columns depending on whether verbose is TRUE or FALSE
#'   respectively.
#' 
#' \describe{
#'
#' \item{\code{subj_id}}{}
#'
#' \item{\code{A}}{Level of within-subject factor A.}
#'
#' \item{\code{B}}{Level of between-subject factor B.}
#'
#' \item{\code{A_c}}{Deviation-coded predictor for A.}
#'
#' \item{\code{B_c}}{Deviation-coded predictor for B.}
#'
#' \item{\code{tnum_r}}{Trial number for the fully randomized version.}
#'
#' \item{\code{tnum_b}}{Trial number for the blocked version.}
#'
#' \item{\code{Y_r}}{Response variable for the fully randomized version.}
#'
#' \item{\code{Y_b}}{Response variable for the fully blocked version.}
#'
#' \item{\code{item_rint}}{Item random intercept effect (verbose mode only.)}
#'
#' \item{\code{item_rslp}}{Item random slope effect (verbose mode only).}
#'
#' \item{\code{subj_rint}}{Subject random intercept effect (verbose mode only.)}
#'
#' \item{\code{subj_rslp}}{Subject random slope effect (verbose mode only).}
#'
#' \item{\code{Y_fit}}{Fitted value (all effects except residual.)}
#' 
#' }
#'
#' @export
sim_2x2 <- function(n_subj = 48, n_items = 48,
                    int = 0, A = 0, B = 0, AB = 0,
                    item_rint = .5, item_rslp = .5, item_rcorr = .5,
                    subj_rint = .5, subj_rslp = .5, subj_rcorr = .5,
                    version = 0L, verbose = FALSE,
                    extra_args = NULL) {

  ## n_obs <- n_items # for future development  
  if ((n_subj %% 4L) || (n_subj <= 0L))
    stop("'n_subj' must be a positive integer that is a multiple of 4")

  if ((n_items %% 2L) || (n_items <= 0L))
    stop("'n_items' must be a positive even integer")

  make_rfx <- function(rint, rslp, rcorr){
    rfx_mx <- matrix(c(rint,
                       sqrt(rint) * sqrt(rslp) * rcorr,
                       sqrt(rint) * sqrt(rslp) * rcorr,
                       rslp), nrow = 2)

    return(rfx_mx)
  }  

  subj_rfx_mx <- make_rfx(subj_rint, subj_rslp, subj_rcorr)
  item_rfx_mx <- make_rfx(item_rint, item_rslp, item_rcorr)
  
  sfx <- as.data.frame(
    MASS::mvrnorm(n_subj,
                  c(subj_rint = 0, subj_rslp = 0),
                  subj_rfx_mx))

  ifx <- as.data.frame(
    MASS::mvrnorm(n_items,
                  c(item_rint = 0, item_rslp = 0),
                  item_rfx_mx))


  sfx[["subj_id"]] <- factor(seq_len(n_subj))
  sfx[["B"]] <- factor(rep(c("B1", "B2"), times = n_subj / 2L))
  sfx[["blk_ord"]] <- factor(rep(rep(1:2, each = 2L), times = n_subj / 4L))

  ifx[["item_id"]] <- factor(seq_len(n_items))
  ifx[["A"]] <- factor(rep(c("A1", "A2"), times = n_items / 2L))

  ## Get complete list of all subject and item combinations
  trials <- expand.grid(
    subj_id = sfx[["subj_id"]],
    item_id = ifx[["item_id"]]
  )

  ## Join subject and item rfx
  trials_subj <- merge(sfx, trials, by = "subj_id")
  dat <- merge(ifx, trials_subj, by = "item_id")

  dat[["A_c"]] <- ifelse(dat[["A"]] == "A1", -.5, .5)
  dat[["B_c"]] <- ifelse(dat[["B"]] == "B1", -.5, .5)

  dat[["Y_fit"]] <-
    int + dat[["subj_rint"]] + dat[["item_rint"]] +
    (A + dat[["subj_rslp"]]) * dat[["A_c"]] +
    (B + dat[["item_rslp"]]) * dat[["B_c"]] +
    AB * dat[["A_c"]] * dat[["B_c"]]

  ds <- split(dat, dat[["subj_id"]])

  ## new vers vvv
  errs <- if (is.numeric(version)) {
    replicate(n_subj, errsim(n_items, version), simplify = FALSE)
  } else {
    rres <- do.call(version, c(list(n_subj = n_subj, n_obs = n_items),
                               extra_args))
    lapply(rres, function(vv) (vv - mean(vv)) / sd(vv))
  }
  
  if (length(errs) != n_subj) {
    stop("user-defined residual function must return a list with length 'n_subj'")
  }

  nobstest <- sapply(errs, length)
  if ( (length(unique(nobstest)) != 1L) || any(unique(nobstest) != n_items) ) {
    stop("all elements in list returned by user-defined residual function must be of length 'n_items'")
  }
  ## new vers ^^^  
  ## errs <- replicate(n_subj, errsim(n_items, version), simplify = FALSE)

  ## add in the errors
  derr <- mapply(function(.d, .e) {
    .d[["tnum_r"]] <- sample(seq_len(nrow(.d)))
    .d2 <- split(.d, .d[["A"]])
    if (.d[["blk_ord"]][1] == 1L) {
      .d2[[1]][["tnum_b"]] <- sample(seq_len(n_items / 2L))
      .d2[[2]][["tnum_b"]] <- sample(seq_len(n_items / 2L)) + (n_items / 2L)
    } else {
      .d2[[1]][["tnum_b"]] <- sample(seq_len(n_items / 2L)) + (n_items / 2L)
      .d2[[2]][["tnum_b"]] <- sample(seq_len(n_items / 2L))
    }
    .d3 <- rbind(.d2[[1]], .d2[[2]])
    .d3[["Y_r"]] <- .d3[["Y_fit"]] + .e[.d3[["tnum_r"]]]
    .d3[["Y_b"]] <- .d3[["Y_fit"]] + .e[.d3[["tnum_b"]]]
    .d3[["err"]] <- .e
    .d3[["tnum_e"]] <- seq_along(.e)
    .d3
  }, ds, errs, SIMPLIFY = FALSE)

  dall <- do.call("rbind", derr)
  rownames(dall) <- NULL

  cols <- c("subj_id", "item_id", "A", "B", "A_c", "B_c", "tnum_r", "tnum_b",
            "Y_r", "Y_b")

  if (verbose) {
    cols <- c(cols, "subj_rint", "subj_rslp", "item_rint", "item_rslp", "Y_fit",
              "tnum_e", "err")
  }

  dall[, cols]
}
