## miscellaneous functions

#' @title Permute residuals
#'
#' Independently permute the residuals for each subject in the data set.
#'
#' @param data The dataset.
#' @param subject Unquoted name of the variable identifying individual
#'   subjects.
#' @param residuals Unquoted variable name for residuals.
#' @return A table with the \code{residuals} column permuted for each subject.
#' @importFrom magrittr %>%
#' @importFrom rlang !!
#' @importFrom rlang :=
#' @export
permute_resids <- function(data, subject, residuals) {
  subj <- rlang::enquo(subject)
  resid <- rlang::enquo(residuals)
  r2 <- rlang::quo_name(resid)
  data %>%
    dplyr::group_by(!!subj) %>%
    dplyr::mutate(!!r2 := sample(!!resid)) %>%
    dplyr::ungroup()
}
