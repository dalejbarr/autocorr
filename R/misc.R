#' Consolidate Simulation Results
#'
#' Consolidate results from multiple simulations into a smaller set of
#' files, removing date/time, host, and PID information.
#'
#' @param inpath Path to directory containing simulation results.
#'
#' @param outpath Path to directory where results files will be stored.
#' 
#' @export
consolidate_results <- function(inpath, outpath) {
  if (!dir.exists(inpath)) {
    stop("directory '", inpath, "' does not exist")
  }
  if (!dir.exists(outpath)) {
    stop("directory '", outpath, "' does not exist; you must create it.\n",
         "  dir.create('", outpath, "')")
  } else {
    if (length(dir(outpath) > 0L))
      stop("directory '", outpath, "' must be empty.")
  }
  allfiles <- dir(inpath, "^ac.*\\.rds$")
  cfgs <- sapply(strsplit(allfiles, "_"), function(.x) {
    case <- if (.x[1] == "ac") {
              as.character(as.integer(.x[12]) - 1L)
            } else {
              .x[12]
            }
    paste(c(.x[3:11], case), collapse = "_")
  })
  files <- unique(cfgs)
  for (i in seq_along(files)) {
    ofile <- file.path(sub("/$", "", outpath),
                       paste0(files[i], ".rds"))
    res <- lapply(file.path(sub("/$", "", inpath),
                            allfiles[cfgs == files[i]]), readRDS)
    res2 <- do.call("rbind", res)
    message("wrote ", sum(cfgs == files[i]), " files (",
            dim(res2)[1], " runs) to ", ofile)
    saveRDS(res2, ofile)
  }
  invisible(NULL)
}
