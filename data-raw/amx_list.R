## This script process the data in sim_results_raw/sine_*/ and
## creates data objects in sim_results/
## 
## sine_phase_results :
##   results from sine-wave simulations with varying phase angle
##
## sine_amp_results :
##   results from sine-wave simulations with varying amplitudes
##
## need to download the directories:
##  https://github.com/dalejbarr/autocorrelation/tree/master/simulation_results/sine_amp
##  https://github.com/dalejbarr/autocorrelation/tree/master/simulation_results/sine_phase
##
## put these two directories as subdirectories in your working directory.
##
library("abind")
library("tidyverse")

process_files <- function(path) {
  all_files <- list.files(path, "^rmx.+\\.rds$", full.names = TRUE)
  finf <- map_int(all_files, ~ as.integer(file.info(.x)$size))
  files_todo <- all_files[finf > 0] # drop empty files

  dat <- tibble(fname = files_todo,
                A = sub(".*rmx_([0-9\\.]+)_.+", "\\1", fname) %>%
                  as.numeric(),
                B = sub(".*rmx_[0-9\\.]+_([0-9\\.]+)_.+", "\\1", fname) %>%
                  as.numeric(),
                AB = sub(".*rmx_[0-9\\.]+_[0-9\\.]+_([0-9\\.]+)_.+", "\\1",
                         fname) %>% as.numeric(),
                nmc = sub(".*rmx_.+_([0-9]{5})_.+", "\\1", fname) %>%
                  as.integer(),
                host = sub(".+_([A-Za-z\\-]+)_[0-9]{3}\\.rds$", "\\1", fname),
                ix = sub(".+_([0-3]{3})\\.rds$", "\\1", fname) %>%
                  as.integer())

  nmc_runs <- dat %>%
    group_by(A, B, AB) %>%
    summarize(tot = sum(nmc)) %>%
    pull(tot) %>%
    unique()

  stopifnot(length(nmc_runs) == 1L)

  fff <- split(dat, paste(dat$A, dat$B, dat$AB))

  map(fff, function(x) {
    amx <- do.call(abind,
                   c(map(x$fname, ~ readRDS(.x)),
                     list(along = 4)))
  })
}

sine_phase_results <- process_files("../sine_phase")
sine_amp_results <- process_files("../sine_amp")

saveRDS(sine_phase_results, file = "../
save(

tools::resaveRdaFiles(f1)
tools::checkRdaFiles(f1)

devtools::use_data(, overwrite = TRUE)
