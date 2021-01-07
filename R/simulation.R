sawtooth <- function(x) {
  n <- sample(2:9, 1)
  ## n = 1 -> sine function
  ## as n increase, becomes more sawtooth like
  ## has 2pi frequency
  sin.k <- function(n,k,x){
    choose(2*n,n-k)/choose(2*n,n)*sin(k*x)/k
  }
  
  rowSums(purrr::map_dfc(setNames(1:n,1:n), ~ sin.k(n,.x,x)))
}

## This function written by R. Harald Baayen
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
#' 
#'     \item{1}{Sine wave with fixed amplitude, varying phase, and white noise.}
#' 
#'     \item{2}{Sine wave with fixed phase, varying amplitude, and white noise.}
#' 
#'     \item{3}{Random walk, high frequency, generated using
#'       [stat_gp()] with gamma = 1 and sigma = 1.}
#'     \item{4}{Random walk, mid frequency, generated using
#'       [stat_gp()] with gamma = 2 and sigma = 1.}
#' 
#'     \item{5}{Multi-scale: combination of 1 and 3.}
#'     \item{6}{Multi-scale: combination of 1 and 4.}
#'     \item{7}{Multi-scale: combination of 2 and 3.}
#'     \item{8}{Multi-scale: combination of 2 and 4.}
#' 
#'     \item{9}{A wiggly function generated from
#'     [mgcv::smooth.construct()], no error}
#' 
#'     \item{10}{A wiggly function generated from
#'     [mgcv::smooth.construct()], with 10 percent white noise}
#'
#'     \item{11}{Analogous to 1, but with a sawtooth pattern instead
#'     of sine.}
#' 
#'     \item{12}{Analogous to 2, but with a sawtooth pattern instead
#'     of sine.}
#' 
#'     \item{13}{Multi-scale: combination of 11 and 3.}
#'     \item{14}{Multi-scale: combination of 11 and 4.}
#'     \item{15}{Multi-scale: combination of 12 and 3.}
#'     \item{16}{Multi-scale: combination of 12 and 4.}
#' }
#' 
#' @return A vector of simulated observations guaranteed to have a
#'   mean of 0 and a standard deviation of 1.
#' 
#' @export
errsim <- function(n_obs, version) {
## errsim <- function(n_obs, version, det_trend, m, m_sw)  {
  
##  if (m_sw==TRUE) {
##    m <- runif(1L,-m,m)
##    }
  
##  det.tr <- function(x,m) {det.trend(det_trend)(x) + m}
  if (n_obs > length(stat_gp(1, 1)$GP))
    stop("'n_obs' must be smaller than ",
         length(stat_gp(1, 1)$GP))  
  
  version_int <- as.integer(version)
  if (is.na(version_int))
    stop("'version' must be an integer")
  maxvers <- 16L
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
  } else if (version == 10L) {
    vv1 <- wiggle(n_obs)
    vv2 <- (vv1 - mean(vv1)) / sd(vv1)
    vv <- vv2 + rnorm(n_obs, sd = sqrt(.1))
  } else if (version == 11L) {
    ## analogous to 1 (fixed amp varying phase)
    phase <- runif(1L, -pi, pi)
    vv <- ((sawtooth(x + phase) / sd(sawtooth(x + phase))) +
           rnorm(length(x), 0, sqrt(.1))) / sqrt(1.1)
    attr(vv, "phase") <- phase
  } else if (version == 12L) {
    ## analogous to 2 (fixed phase varying amp)
    ampvar <- runif(1L, .2, .8)
    ## ampvar <- 1
    noisevar <- 1 - ampvar
    ##noisevar <- 0
    vv <- ((sawtooth(x) / sd(sawtooth(x))) * sqrt(ampvar)) +
      rnorm(n_obs, sd = sqrt(noisevar))
    attr(vv, "amp") <- ampvar
  } else if (version == 13L) {
    ## mixed 11 and 3
    det.tramp <- .9
    noise_lvl <- 1 - det.tramp
    phase <- runif(1L, -pi, pi)
    sdat <- sawtooth(x + phase) / sd(sawtooth(x + phase))
    ndat0 <- stat_gp(1, 1)$GP[seq_len(n_obs)]
    ndat1 <- (ndat0 - mean(ndat0))
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(det.tramp) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "phase") <- phase
  } else if (version == 14L) {
    ## mixed 11 and 4
    det.tramp <- .9
    noise_lvl <- 1 - det.tramp
    phase <- runif(1L, -pi, pi)
    sdat <- sawtooth(x + phase) / sd(sawtooth(x + phase))
    ndat0 <- stat_gp(1, 2)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(det.tramp) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "phase") <- phase
  } else if (version == 15L) {
    ## mixed 12 and 3
    ampvar <- runif(1L, .2, .8)
    noise_lvl <- 1 - ampvar
    sdat <- sawtooth(x) / sd(sawtooth(x))
    ndat0 <- stat_gp(1, 1)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(ampvar) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "amp") <- ampvar
  } else if (version == 16L) {
    ## mixed 12 and 4
    ampvar <- runif(1L, .2, .8)
    noise_lvl <- 1 - ampvar
    sdat <- sawtooth(x) / sd(sawtooth(x))
    ndat0 <- stat_gp(1, 2)$GP[seq_len(n_obs)]
    ndat1 <- ndat0 - mean(ndat0)
    ndat <- ndat1 / sd(ndat1)
    vv <- sqrt(ampvar) * sdat +
      sqrt(noise_lvl) * ndat
    attr(vv, "amp") <- ampvar        
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
#' @param version How to generate residuals: either an integer
#'   representing the Scenario number (see \code{\link{errsim}}) or
#'   the name of a user-defined function.
#'
#' @param verbose Whether the data frame should include GLM
#'   components.
#'
#' @param extra_args A list of extra arguments to be passed to a
#'   user-defined function to generate residuals.
#'
#' @details When used with a user-defined function to generate
#'   residuals, the residuals for each subject will be standardized
#'   (i.e., converted to z-scores) before they are combined with other
#'   model components.
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
                    verbose = FALSE,
                    extra_args = NULL) {

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

  errs <- if (is.numeric(version)) {
    replicate(n_subj, errsim(n_obs, version), simplify = FALSE)
  } else {
    rres <- do.call(version, c(list(n_subj = n_subj, n_obs = n_obs), extra_args))
    lapply(rres, function(vv) (vv - mean(vv)) / sd(vv))
  }
  
  if (length(errs) != n_subj) {
    stop("user-defined residual function must return a list with length 'n_subj'")
  }

  nobstest <- sapply(errs, length)
  if ( (length(unique(nobstest)) != 1L) || any(unique(nobstest) != n_obs) ) {
    stop("all elements in list returned by user-defined residual function must be of length 'n_obs'")
  }

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
#' @param bam_args Any other arguments to be passed on to
#'   [mgcv::bam()] (for the fitting of the GAMM models only.)
#'
#' @param fit_blocked Whether to fit the blocked version in addition
#'   to the randomized version.
#'
#' @param fit_lmem Whether to fit the LMEM in addition to the GAMM.
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
                    dontfit = FALSE, m = NA, k = -1, bam_args = NULL,
                    fit_blocked = TRUE, fit_lmem = TRUE) {
  ## function to extract model statistics
  mod_stats <- function(m_gamm, m_lmem = NULL) {
    if (!is.null(m_lmem)) {
      mc <- anova(m_gamm, m_lmem)
      rdf_g <- mc[["Resid. Df"]][1]
      rdev_g <- mc[["Resid. Dev"]][1]
      rdf_l <- mc[["Resid. Df"]][2]
      rdev_l <- mc[["Resid. Dev"]][2]
      cf_l <- coef(m_lmem)[1:4]
      se_l <- sqrt(diag(vcov(m_lmem)[1:4, 1:4]))
      m_lmem_s <- summary(m_lmem)
      p_l <- m_lmem_s[["p.table"]][, "Pr(>|t|)"]
      aic_l <- AIC(m_lmem)
    } else {
      rdf_g <- rdev_g <- rdf_l <- rdev_l <- aic_l <- NA_real_
      cf_l <- se_l <- p_l <- rep(NA_real_, 4)
    }
    m_gamm_s <- summary(m_gamm)
    v <- c(coef(m_gamm)[1:4],
           sqrt(diag(vcov(m_gamm)[1:4, 1:4])),
           m_gamm_s[["p.table"]][, "Pr(>|t|)"],
           AIC(m_gamm),
           rdf_g,
           rdev_g,
           cf_l, ## lmem stats start here
           se_l,
           p_l,
           aic_l,
           rdf_l,
           rdev_l)
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
  f_rand2 <- as.formula(paste0("Y_r ~", form_rhs_no_gamm)) ## lmem
  f_block <- as.formula(paste0("Y_b ~", form_rhs_b))
  f_block2 <- as.formula(paste0("Y_b ~", form_rhs_no_gamm)) ## lmem

  if (dontfit) {
    list(randomized = list(GAM = f_rand, LMM = f_rand2),
         blocked = list(GAM = f_block, LMM = f_block2))
  } else {

    dat_r <- dat[order(dat$subj_id, dat$tnum_r), ]
    dat_r[["first"]] <- dat_r[["tnum_r"]] == 1L
    dat_b <- dat[order(dat$subj_id, dat$tnum_b), ]
    dat_b[["first"]] <- dat_b[["tnum_b"]] == 1L
      
    ## fit the GAMM models using mgcv::bam
    mod_rand <- do.call(getExportedValue("mgcv", "bam"),
                        args = c(list(formula = f_rand,
                                      data = dat_r,
                                      AR.start = dat_r$first), bam_args))

    ms_rand <- mod_stats(mod_rand, NULL)
    
    if (fit_blocked) {
      mod_block <- do.call(getExportedValue("mgcv", "bam"),
                           args = c(list(formula = f_block,
                                         data = dat_b,
                                         AR.start = dat_b$first), bam_args))
      ms_block <- mod_stats(mod_block, NULL)
    } else {
      ms_block <- rep(NA_real_, length(ms_rand))
      names(ms_block) <- names(ms_rand)
    }

    if (fit_lmem) {
      ## fit the LMEM models using mgcv::bam
      mod_rand_2 <- do.call(getExportedValue("mgcv", "bam"),
                            args = list(formula = f_rand2,
                                        data = dat_r))

      ms_rand <- mod_stats(mod_rand, mod_rand_2)
    
      if (fit_blocked) {
        mod_block_2 <- do.call(getExportedValue("mgcv", "bam"),
                               args = list(formula = f_block2,
                                           data = dat_b))
        
        ms_block <- mod_stats(mod_block, mod_block_2)
      }
    }

    array(c(ms_rand,
            ms_block),
          dim = c(15, 2, 2),
          dimnames = list(parm = names(mod_stats(mod_rand, NULL))[1:15],
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
#' @param version Autocorrelation case number; either an integer
#'   corresponding to Scenario number (see \link{errsim}) or the name
#'   of a user-defined function specified as a character string.
#'
#' @param os_always Whether to always fit an overall smooth, or only
#'   for [errsim()] versions 2, 7, 8, 12, 15, and 16. For user-defined
#'   functions, default is to fit a common smooth.
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
#' @param fit_blocked Whether to fit the blocked version in addition
#'   to the randomized version.
#'
#' @param fit_lmem Whether to fit the LMEM in addition to the GAMM.
#'
#' @param outfile Name of output file.
#'
#' @param extra_args Extra args to be passed along to any user-defined
#'   function for generating residuals.
#'
#' @details The behavior of \code{os_always} depends on whether
#'   \code{version} has been specified as a scenario number, in which
#'   case it overrides the scenario-dependent selection of an overall
#'   smooth, always fitting an overall smooth. If version is the name
#'   of a user-defined function, then \code{os_always} determines
#'   whether the model includes an overall smooth.
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
                  os_always = is.character(version),
                  m = NA,
                  k = -1,
                  bam_args = NULL,
                  fit_blocked = TRUE,
                  fit_lmem = TRUE,
                  outfile = sprintf(
                    "acs_%05d_%03d_%03d_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%0.2f_%s_%s_%s_%s.rds",
                    nmc, n_subj, n_obs,
                    A, B, AB,
                    rint_range[1], rint_range[2],
                    rslp_range[1], rslp_range[2],
                    if (is.numeric(version)) sprintf("%02d", version) else version,
                    format(Sys.time(), "%Y-%m-%d-%H-%M-%S"),
                    Sys.info()[["nodename"]],
                    Sys.getpid()),
                  extra_args = NULL) {
  
  tfile <- tempfile(fileext = ".csv")
                    
  cs <- if (os_always) {
          TRUE
        } else {
          if (is.numeric(version)) {
            version %in% c(2L, 7L, 8L, 12L, 15L, 16L)
          } else {
            os_always
          }
        }

  append <- FALSE

  for (i in seq_len(nmc)) {
    rint <- runif(1, rint_range[1], rint_range[2])
    rslp <- runif(1, rslp_range[1], rslp_range[2])
    rcorr <- runif(1, rcorr_range[1], rcorr_range[2])

    dat <- sim_2x2(n_subj = n_subj, n_obs = n_obs,
                   int = 0, A = A, B = B, AB = AB,
                   rint = rint, rslp = rslp, rcorr = rcorr,
                   version = version, extra_args = extra_args)
    res <- fit_2x2(dat = dat, cs = cs, m = m, k = k,
                   bam_args = bam_args, fit_blocked = fit_blocked,
                   fit_lmem = fit_lmem)
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
