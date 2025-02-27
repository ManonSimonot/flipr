#' Linar Regression Permutation Test
#'
#' This function carries out an hypothesis test in which the null hypothesis is
#' that predictors don't have an effect on the response variable against the
#' alternative hypothesis that they do have an effect on the response variable.
#'
#' @section User-supplied statistic function:
#' A user-specified function should have at least two arguments:
#'
#' - the first argument should be a numeric vector representing the response
#' variable.
#' - the second argument should be a data drame of the `n` pooled data points.
#' - the third argument is optional and represent the index of the coefficient
#' to be tested in the case of a test on a specific coefficient.
#'
#' See the [`stat_lm()`] function for an example.
#'
#' @param original_response An numeric vector representing the response variable.
#' @param fitted The fitted values of the reduced model.
#' @param residuals The residual values of the reduced model.
#' @param vars A data frame of all variables used as predictors in the model.
#' @param index For tests on one coefficient, an integer specifying the index of
#' the coefficient to be tested.
#' @param stats A list of functions produced by [`rlang::as_function()`]
#'   specifying the chosen test statistic(s). A number of test statistic
#'   functions are implemented in the package and can be used as such.
#'   Alternatively, one can provide its own implementation of test statistics
#'   that (s)he deems relevant for the problem at hand. See the section
#'   *User-supplied statistic function* for more information on how these
#'   user-supplied functions should be structured for compatibility with the
#'   **flipr** framework. Defaults to `list(stat_lm)`.
#' @param B The number of sampled permutations. Default is `1000L`.
#' @param M The total number of possible permutations. Defaults to `NULL`, which
#'   means that it is automatically computed from the given sample size(s).
#' @param alternative A single string or a character vector specifying whether
#'   the p-value is right-tailed, left-tailed or two-tailed. Choices are
#'   `"right_tail"`, `"left_tail"` and `"two_tail"`. Default is `"two_tail"`. If
#'   a single string is provided, it is assumed that it should be applied to all
#'   test statistics provided by the user. Alternative, the length of
#'   `alternative` should match the length of the `stats` parameter and it is
#'   assumed that there is a one-to-one correspondence. Defaults to
#'   `"right_tail"`.
#' @param combine_with A string specifying the combining function to be used to
#'   compute the single test statistic value from the set of p-value estimates
#'   obtained during the non-parametric combination testing procedure. For now,
#'   choices are either `"tippett"` or `"fisher"`. Defaults to `"tippett"`,
#'   which picks Tippett's function.
#' @param type A string specifying which formula should be used to compute the
#'   p-value. Choices are `exact`, `upper_bound` and `estimate`. See Phipson &
#'   Smith (2010) for details. Defaults to `"exact"`.
#' @param seed An integer specifying the seed of the random generator useful for
#'   result reproducibility or method comparisons. Default is `NULL`.
#' @param ... Extra parameters specific to some statistics.
#'
#' @return A [`base::list`] with 4 components:
#' - `observed`: the value of the (possible combined) test statistic(s) using
#' the original memberhips of data points;
#' - `pvalue`: the permutation p-value;
#' - `null_distribution`: the values of the (possible combined) test statistic(s)
#' using the permuted memberhips of data points;
#' - `permutations`: the permutations that were effectively sampled to produce
#' the null distribution.
#' @export
#'
#' @examples
#' Y <- iris$Sepal.Length
#' regressors <- iris[2:5]
#' reduced_vars <- regressors[- c(1)]
#' fit <- lm(Y ~ ., data = reduced_vars)
#' out2 <- regression_test(
#'   fitted <- fit$fitted,
#'   residuals <- fit$residuals,
#'   vars = iris[2:5],
#'   index = 1,
#'   stats = list(stat_lm)
#' )
#' out2$pvalue
regression_test <- function(fitted, residuals,
                            vars, indexes,
                            stats = list(stat_lm),
                            B = 1000L,
                            M = NULL,
                            alternative = "right_tail",
                            combine_with = "tippett",
                            type = "exact",
                            seed = NULL,
                            ...) {

  if (rlang::is_bare_numeric(vars) || is.matrix(vars) || is.data.frame(vars))
    vars <- list(vars)

  if (!is.list(vars))
    abort("The {.arg data} argument should be either of class {.cls numeric} or
          of class {.cls matrix} or of class {.cls list}}.")

  response <- fitted + residuals
  n <- length(response)

  # Compute total number of permutations yielding to distinct
  # values of the test statistic
  if (is.null(M))
    M <- factorial(n) - 1

  # Case of global test
  if (indexes == 0) {
    # Generate permutation data by permuting response
    perm_data <- replicate(B, sample(response, size = n))
    perm_data <- cbind(response, perm_data)
  }

  # Case of one coefficient test
  else {
    # Generate permutation data by permuting residuals of reduced model
    perm_residuals <- replicate(B, sample(residuals, size = n))

    # Permuted response
    perm_data <- fitted + perm_residuals
    perm_data <- cbind(response, perm_data)
  }

  run_permutation_scheme(
    type = type,
    alternative = alternative,
    stats = stats,
    B = B,
    perm_data = perm_data,
    stat_data = vars,
    indexes = indexes,
    M = M,
    combine_with = combine_with,
    ...
  )
}
