#' Test Statistics for the Linear Regression Problem
#'
#' This is a collection of functions that provide test statistics to be used
#' into the permutation scheme for performing Linear Regression. These test
#' statistics can be divided into two categories: traditional statistics that
#' use empirical moments and inter-point statistics that only rely on pairwise
#' dissimilarities between data points.
#'
#' @section Traditional Test Statistics:
#'
#' - [`stat_lm_global()`] implements the F statistic used to test
#' globally the coefficients of the regression.
#'
#' - [`stat_lm()`] implements the t statistic used to test one
#' coefficient of the regression.
#'
#' @param vars A data frame of all variables used as predictors in the model.
#' @param response An numeric vector representing the response variable.
#' @param index For tests on one coefficient, an integer specifying the index of
#' the coefficient to be tested.
#' @param ... Extra parameters specific to some statistics.
#'
#' @return A numeric value storing the value of test statistic given the
#'   (possibly permuted) response specified by `response`.
#' @name regression-stats
#'
#' @examples
#' response_var <- iris$Sepal.Length
#' predictors <- iris[2:5]
#' stat_lm_global(predictors, response_var)
#' stat_lm(predictors, response_var, index = 1)
#'
NULL

#' @rdname regression-stats
#' @export
stat_lm_global <- function(vars, response, ...) {
  response_var <- response
  vars <- as.data.frame(vars)
  fit <- stats::lm(response_var ~ ., data = vars)
  summary(fit)$f[1] #F-statistic
}

#' @rdname regression-stats
#' @export
stat_lm <- function(vars, response, index, ...) {
  response_var <- response
  if (!is.data.frame(vars)) { vars <- as.data.frame(vars) }
  fit <- stats::lm(response_var ~ ., data = vars)
  abs(summary(fit)$coefficients[index + 1, 3]) # t-statistic
}


