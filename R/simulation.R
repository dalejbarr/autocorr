#' This function written by R. Harald Baayen
wiggle <- function(max_time, stddev=2, k=10, scaling=1) {
  x = seq(0, 1, length = max_time)
                                        # time points
  sm <- mgcv::smooth.construct(mgcv::s(x, k = k, bs = "tp"),
                               data = data.frame(x = x),
                               knots = NULL)

  b = rnorm(k, 0, stddev) * scaling  # simulate coefficients sampled from
                                     # an N(0, stdev) distribution;
  X <- b*t(sm$X)                     # weight the basis functions
  f<-X[1,]*0
                                        # and add them up
  for (i in 1:k) {
    f <- f + X[i,]
  }
  return(f)
}


#' Simulate Autocorrelated Errors
#'
#' @param n_obs Number of simulated observations. Must be less than
#'   275 (length of output from [stat_gp()] with `gamma = 1`, `sigma = 1`.)
#' 
#' @param version An integer specifying the error autocorrelation
#'   case number; one of the following.
#' 
#'   \describe{
#'     \item{0}{No autocorrelation (white noise).}
#'     \item{1}{Sine wave with fixed amplitude, varying phase, and white noise.}
#'     \item{2}{Sine wave with fixed phase, varying amplitude, and white noise.}
#'     \item{3}{Random walk, high frequency, generated using
#'       [stat_gp()] with gamma = 1 and sigma = 1.}
#'     \item{4}{Random walk, mid frequency, generated using
#'       [stat_gp()] with gamma = 2 and sigma = 1.}
#'     \item{5}{Multi-scale: combination of 1 and 3.}
#'     \item{6}{Multi-scale: combination of 1 and 4.}
#'     \item{7}{Multi-scale: combination of 2 and 3.}
#'     \item{8}{Multi-scale: combination of 2 and 4.}
#'     \item{9}{GAMM: a wiggly function generated from a GAMM}
#' }
#' 
#' @return A vector of simulated observations guaranteed to have a
#'   mean of 0 and a standard deviation of 1.
#' 
#' @export
errsim <- function(n_obs, version) {

  if (n_obs > length(stat_gp(1, 1)$GP))
    stop("'n_obs' must be smaller than ",
         length(stat_gp(1, 1)$GP))
  
  version_int <- as.integer(version)
  if (is.na(version_int))
    stop("'version' must be an integer")
  maxvers <- 9L
  if ((version_int < 0L) || (version_int > maxvers))
    stop("'version' must be between 0 and ", maxvers)
  
  x <- seq(-pi, pi, length.out = n_obs)

  if (version == 0L) {
    ## white noise
    vv <- rnorm(n_obs)
  } else if (version == 1L) {
    ## varying phase + white noise
    phase <- runif(1L, -pi, pi)
    vv <- ((sin(x + phase) / sd(sin(x + phase))) +
           rnorm(length(x), 0, sqrt(.1))) / sqrt(1.1)
    attr(vv, "phase") <- phase
  } else if (version == 2L) {
    ## fixed phase, varying amp + white noise
    ampvar <- runif(1L, .2, .8)
    noisevar <- 1 - ampvar
    vv <- ((sin(x) / sd(sin(x))) * sqrt(ampvar)) +
      rnorm(n_obs, sd = sqrt(noisevar))
    attr(vv, "amp") <- ampvar
  } else if (version == 3L) {
    ## random walk, high frequency
    vv <- stat_gp(1, 1)$GP[seq_len(n_obs)]
  } else if (version == 4L) {
    ## random walk, mid frequency
    vv <- stat_gp(1, 2)$GP[seq_len(n_obs)]
  } else if (version == 5L) {
    ## phase plus high freq random walk
    sinamp <- .9
    noise_lvl <- 1 - sinamp
    phase <- runif(1L, -pi, pi)
    sdat <- sin(x + phase) / sd(sin(x + phase))
    ndat0 <- stat_gp(1, 1)$GP[seq_len(n_obs)]
    ndat1 <- (ndat0 - mean(ndat0))
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(sinamp) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "phase") <- phase
  } else if (version == 6L) {
    ## phase plus mid freq random walk
    sinamp <- .9
    noise_lvl <- 1 - sinamp
    phase <- runif(1L, -pi, pi)
    sdat <- sin(x + phase) / sd(sin(x + phase))
    ndat0 <- stat_gp(1, 2)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(sinamp) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "phase") <- phase
  } else if (version == 7L) {
    ## varying amp plus high freq random walk
    ampvar <- runif(1L, .2, .8)
    noise_lvl <- 1 - ampvar
    sdat <- sin(x) / sd(sin(x))
    ndat0 <- stat_gp(1, 1)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(ampvar) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "amp") <- ampvar
  } else if (version == 8L) {
    ## varying amp plus mid freq random walk
    ampvar <- runif(1L, .2, .8)
    noise_lvl <- 1 - ampvar
    sdat <- sin(x) / sd(sin(x))
    ndat0 <- stat_gp(1, 2)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(ampvar) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "amp") <- ampvar    
  } else if (version == 9L) {
    vv <- wiggle(n_obs)
  } else {
    stop("version '", version, "' not recognized")
  }

  (vv - mean(vv)) / sd(vv)
}

#' Simulate 2x2 data with autocorrelated errors
#'
#' @param n_subj Number of subjects to simulate. Must be a positive
#'   integer that is a multiple of 2.
#' 
#' @param n_obs Number of observations per subject. Must be a positive
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
#' @param rint Random intercept variance.
#'
#' @param rslp Random slope variance (for factor A).
#'
#' @param rcorr Random correlation.
#'
#' @param version Autocorrelation case number; an integer between 0
#'   and 8 (see \code{\link{errsim}}).
#'
#' @param verbose Whether the data frame should GLM components.
#' 
#' @return A data frame with \code{n_subj * n_obs} rows and either 9
#'   or 12 columns depending on whether verbose is TRUE or FALSE
#'   respectively.
#' 
#'   \describe{
#'     \item{\code{subj_id}}{}
#'     \item{\code{A}}{Level of within-subject factor A.}
#'     \item{\code{B}}{Level of between-subject factor B.}
#'     \item{\code{A_c}}{Deviation-coded predictor for A.}
#'     \item{\code{B_c}}{Deviation-coded predictor for B.}
#'     \item{\code{tnum_r}}{Trial number for the fully randomized version.}
#'     \item{\code{tnum_b}}{Trial number for the blocked version.}
#'     \item{\code{Y_r}}{Response variable for the fully randomized version.}
#'     \item{\code{Y_b}}{Response variable for the fully blocked version.}
#'     \item{\code{rint}}{Random intercept effect (verbose mode only.)}
#'     \item{\code{rslp}}{Random slope effect (verbose mode only).}
#'     \item{\code{Y_fit}}{Fitted value (all effects except residual.)}
#' }
#'
#' @export
sim_2x2 <- function(n_subj = 48, n_obs = 48,
                    int = 0, A = 0, B = 0, AB = 0,
                    rint = .5, rslp = .5, rcorr = .5,
                    version = 0L,
                    verbose = FALSE) {

  if ((n_subj %% 4L) || (n_subj <= 0L))
    stop("'n_subj' must be a positive integer that is a multiple of 4")

  if ((n_obs %% 2L) || (n_obs <= 0L))
    stop("'n_obs' must be a positive even integer")
  
  rfx_mx <- matrix(c(rint,
                     sqrt(rint) * sqrt(rslp) * rcorr,
                     sqrt(rint) * sqrt(rslp) * rcorr,
                     rslp), nrow = 2)

  sfx <- as.data.frame(
    MASS::mvrnorm(n_subj, c(rint = 0, rslp = 0), rfx_mx))
  sfx[["subj_id"]] <- factor(seq_len(n_subj))
  sfx[["B"]] <- factor(rep(c("B1", "B2"), times = n_subj / 2L))
  sfx[["blk_ord"]] <- factor(rep(rep(1:2, each = 2L), times = n_subj / 4L))

  trials <- expand.grid(
    subj_id = factor(seq_len(n_subj)),
    A = factor(rep(c("A1", "A2"), each = n_obs / 2L)))

  dat <- merge(sfx, trials, by = "subj_id")
  dat[["A_c"]] <- ifelse(dat[["A"]] == "A1", -.5, .5)
  dat[["B_c"]] <- ifelse(dat[["B"]] == "B1", -.5, .5)
  dat[["Y_fit"]] <-
    int + dat[["rint"]] +
    (A + dat[["rslp"]]) * dat[["A_c"]] +
    B * dat[["B_c"]] +
    AB * dat[["A_c"]] * dat[["B_c"]]

  ds <- split(dat, dat[["subj_id"]])
  
  errs <- replicate(n_subj, errsim(n_obs, version), simplify = FALSE)

  ## add in the errors
  derr <- mapply(function(.d, .e) {
    .d[["tnum_r"]] <- sample(seq_len(nrow(.d)))
    .d2 <- split(.d, .d[["A"]])
    if (.d[["blk_ord"]][1] == 1L) {
      .d2[[1]][["tnum_b"]] <- sample(seq_len(n_obs / 2L))
      .d2[[2]][["tnum_b"]] <- sample(seq_len(n_obs / 2L)) + (n_obs / 2L)
    } else {
      .d2[[1]][["tnum_b"]] <- sample(seq_len(n_obs / 2L)) + (n_obs / 2L)
      .d2[[2]][["tnum_b"]] <- sample(seq_len(n_obs / 2L))
    }
    .d3 <- rbind(.d2[[1]], .d2[[2]])
    .d3[["Y_r"]] <- .d3[["Y_fit"]] + .e[.d3[["tnum_r"]]]
    .d3[["Y_b"]] <- .d3[["Y_fit"]] + .e[.d3[["tnum_b"]]]
    .d3
  }, ds, errs, SIMPLIFY = FALSE)

  dall <- do.call("rbind", derr)
  rownames(dall) <- NULL

  cols <- c("subj_id", "A", "B", "A_c", "B_c", "tnum_r", "tnum_b", "Y_r", "Y_b")

  if (verbose) {
    cols <- c(cols, "rint", "rslp", "Y_fit")
  }

  dall[, cols]
}

#' Fit models to 2x2 data with autocorrelated errors
#'
#' @param dat Data generated using [sim_2x2()].
#' 
#' @param cs Fit smooths as main effect of trial (in addition to
#'   by-subject factor smooths) for the GAMM models?
#' 
#' @param by_subj_fs Include by-subject factor smooths in GAMM?
#'
#' @param dontfit If \code{TRUE}, don't fit the models, just return
#'   the model formulas. Used for debugging.
#'
#' @param m The `m` parameter to be passed on to any factor smooths
#'   (specified in the [mgcv::s()] function).
#'
#' @param k The `k` parameter to be passed on to any factor smooths
#'   (specified in the [mgcv::s()] function).
#'
#' @param bam_args Any other arguments to be passed on to [mgcv::bam()].
#'
#' @details Fits four models using \link[mgcv]{bam}:
#'
#' \describe{
#'   \item{1}{A Generalized Additive Mixed Model (GAMM) for the blocked DV (`Y_b`);}
#'   \item{2}{A Linear Mixed-Effects Model (LMEM) for the blocked DV (`Y_b`);}
#'   \item{3}{A GAMM for the randomized DV (`Y_r`); and}
#'   \item{4}{A LMEM for the randomized DV (`Y_r`).}
#' }
#' 
#' @return 15x2x2 array with model statistics (or just the model
#'   formulas if `dontfit` is `TRUE`). Second dimension is
#'   GAMM or LMEM, third is randomized or blocked, and first dimension
#'   is as follows.
#'
#' \describe{
#' 
#'   \item{`e_int`}{Estimated intercept.}
#'
#'   \item{`e_A`}{Estimated main effect of A.}
#'
#'   \item{`e_B`}{Estimated main effect of B.}
#'
#'   \item{`e_AB`}{Estimated AB interaction.}
#'
#'   \item{`se_int`}{Standard error for the intercept.}
#'
#'   \item{`se_A`}{Standard error for A.}
#'
#'   \item{`se_B`}{Standard error for B.}
#'
#'   \item{`se_AB`}{Standard error for AB.}
#'
#'   \item{`p_int`}{P-value for the intercept.}
#'
#'   \item{`p_A`}{P-value for the main effect of A.}
#'
#'   \item{`p_B`}{P-value for the main effect of B.}
#' 
#'   \item{`p_AB`}{P-value for the AB interaction.}
#'
#'   \item{`AIC`}{Akaike Information Criterion for the model.}
#'
#'   \item{`resid_df`}{Residual degrees of freedom.}
#'
#'   \item{`resid_dev`}{Residual deviance.}
#' }
#' 
#' @export
fit_2x2 <- function(dat, cs = FALSE, by_subj_fs = TRUE,
                    dontfit = FALSE, m = NA, k = -1, bam_args = NULL) {
  ## function to extract model statistics
  mod_stats <- function(m_gamm, m_lmem) {
    mc <- anova(m_gamm, m_lmem)
    m_gamm_s <- summary(m_gamm)
    m_lmem_s <- summary(m_lmem)
    v <- c(coef(m_gamm)[1:4],
           sqrt(diag(vcov(m_gamm)[1:4, 1:4])),
           m_gamm_s[["p.table"]][, "Pr(>|t|)"],
           AIC(m_gamm), mc[["Resid. Df"]][1],
           mc[["Resid. Dev"]][1],
           coef(m_lmem)[1:4],
           sqrt(diag(vcov(m_lmem)[1:4, 1:4])),
           m_lmem_s[["p.table"]][, "Pr(>|t|)"],
           AIC(m_lmem), mc[["Resid. Df"]][2],
           mc[["Resid. Dev"]][2])
    vn <- c("e_int", "e_A", "e_B", "e_AB",
            "se_int", "se_A", "se_B", "se_AB",
            "p_int", "p_A", "p_B", "p_AB",
            "AIC", "resid_df", "resid_dev")
    names(v) <- rep(vn, 2)
    v
  }
  ## determine the model formula
  form_rhs <- "A_c * B_c + \n   s(subj_id, A_c, bs = \"re\")"
  form_rhs_b <- form_rhs # blocked version
  
  form_rhs_no_gamm <- paste0(form_rhs,
                             " + \n   s(subj_id, bs = \"re\")")
  
  if (cs) {
    form_rhs <- paste0(form_rhs,
                       " + \n   s(tnum_r, bs = \"tp\")")
    form_rhs_b <- paste0(form_rhs_b,
                         " + \n   s(tnum_b, bs = \"tp\")")
  }

  if (by_subj_fs) {
    form_rhs <- paste(form_rhs,
                      sprintf(" + \n   s(tnum_r, subj_id, k = %d, m = %d, bs = \"fs\")",
                              k, m))
    form_rhs_b <- paste(form_rhs_b,
                        sprintf(" + \n   s(tnum_b, subj_id, k = %d, m = %d, bs = \"fs\")",
                                k, m))
  }

  f_rand <- as.formula(paste0("Y_r ~", form_rhs))
  f_rand2 <- as.formula(paste0("Y_r ~", form_rhs_no_gamm))
  f_block <- as.formula(paste0("Y_b ~", form_rhs_b))
  f_block2 <- as.formula(paste0("Y_b ~", form_rhs_no_gamm))

  if (dontfit) {
    list(randomized = list(GAM = f_rand, LMM = f_rand2),
         blocked = list(GAM = f_block, LMM = f_block2))
  } else {  
    ## fit the GAMM models using mgcv::bam
    mod_rand <- do.call(getExportedValue("mgcv", "bam"),
                        args = c(list(formula = f_rand,
                                      data = dat), bam_args))
    mod_block <- do.call(getExportedValue("mgcv", "bam"),
                         args = c(list(formula = f_block,
                                       data = dat), bam_args))

    ## fit the LMEM models using mgcv::bam
    mod_rand_2 <- do.call(getExportedValue("mgcv", "bam"),
                          args = c(list(formula = f_rand2,
                                        data = dat), bam_args))
    mod_block_2 <- do.call(getExportedValue("mgcv", "bam"),
                          args = c(list(formula = f_block2,
                                        data = dat), bam_args))

    array(c(mod_stats(mod_rand, mod_rand_2),
            mod_stats(mod_block, mod_block_2)),
          dim = c(15, 2, 2),
          dimnames = list(parm = names(mod_stats(mod_rand, mod_rand_2))[1:15],
                          mod = c("GAMM", "LMEM"),
                          vers = c("randomized", "blocked")))
  }
}

#' Monte Carlo Simulation of Data Analysis with Trial-by-Trial Variation
#'
#' @param nmc Number of Monte Carlo runs.
#' 
#' @param n_subj Number of subjects to simulate. Must be a positive
#'   integer that is a multiple of 2.
#' 
#' @param n_obs Number of observations per subject. Must be a positive
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
#' @param rint_range Range of random intercept variance.
#'
#' @param rslp_range Range of random slope variance (for factor A).
#'
#' @param rcorr_range Range of random correlation.
#'
#' @param version Autocorrelation case number; an integer between 0 and 8 (see \link{errsim}).
#'
#' @param os_always Whether to always fit an overall smooth, or only
#'   for [errsim()] versions 2, 7, and 8.
#' 
#' @param m The `m` parameter to be passed on to any factor smooths
#'   (specified in the [mgcv::s()] function).
#'
#' @param k The `k` parameter to be passed on to any factor smooths
#'   (specified in the [mgcv::s()] function).
#'
#' @param bam_args List of additional arguments to be passed onto
#'   [mgcv::bam()].
#' 
#' @param outfile Name of output file.
#'
#' @return Returns NULL.
#' 
#' @export
mcsim <- function(nmc,
                  n_subj = 48L,
                  n_obs = 48L,
                  A = 0, B = 0, AB = 0,
                  rint_range = blst_quantiles()[, "subj_int"],
                  rslp_range = blst_quantiles()[, "subj_slp"],
                  rcorr_range = c(-.8, .8),
                  version = 0L,
                  os_always = FALSE,
                  m = NA,
                  k = -1,
                  bam_args = NULL,
                  outfile = sprintf(
                    "acs_%05d_%03d_%03d_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%02d_%s_%s_%s.rds",
                    nmc, n_subj, n_obs,
                    A, B, AB,
                    rint_range[1], rint_range[2],
                    rslp_range[1], rslp_range[2],
                    version,
                    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
                    Sys.info()[["nodename"]],
                    Sys.getpid())) {
  
  tfile <- tempfile(fileext = ".csv")
                    
  cs <- if (os_always) {
          TRUE
        } else {
          version %in% c(2L, 7L, 8L)
        }
  
  append <- FALSE

  for (i in seq_len(nmc)) {
    rint <- runif(1, rint_range[1], rint_range[2])
    rslp <- runif(1, rslp_range[1], rslp_range[2])
    rcorr <- runif(1, rcorr_range[1], rcorr_range[2])

    dat <- sim_2x2(n_subj = n_subj, n_obs = n_obs,
                   int = 0, A = A, B = B, AB = AB,
                   rint = rint, rslp = rslp, rcorr = rcorr,
                   version = version)
    res <- fit_2x2(dat = dat, cs = cs, m = m, k = k,
                   bam_args = bam_args)
    vv <- c(rint, rslp, rcorr, res)
    readr::write_csv(as.data.frame(as.list(vv)),
                     tfile,
                     append = TRUE, col_names = FALSE)
  }

  cnames <- paste0(rep(c("G.r.", "L.r.", "G.b.", "L.b."),
                       each = length(dimnames(res)[[1]])),
                   dimnames(res)[[1]])
  
  ff <- readr::read_csv(tfile, c("rint", "rslp", "rcorr", cnames),
                        col_types = readr::cols(.default = readr::col_double()))
  file.remove(tfile)

  if (!is.null(outfile)) {
    saveRDS(ff, outfile)
    message(outfile)
  } else {
    invisible(ff)
  }
}
