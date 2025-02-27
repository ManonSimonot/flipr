compute_confidence_interval <- function(pf,
                                        point_estimate,
                                        conf_level = 0.95,
                                        type = c("interval", "lower_bound", "upper_bound"),
                                        regression = FALSE) {
  if (pf$nparams != 1)
    abort("This function computes a confidence interval. It therefore expect a one-dimensional plausibility function.")

  type <- match.arg(type)

  if (regression) {
    cost <- function(.x) {
      key <- names(pf$parameters)
      named_params <- list()
      named_params[[key]] <- .x
      pvalues <- pf$get_value(named_params)
      pvalues - (1 - conf_level)
    }
  } else {
    cost <- function(.x) {
      pvalues <- pf$get_value(.x)
      pvalues - (1 - conf_level)
    }
  }

  if (point_estimate > 0) {
    lbi <- c(point_estimate / 10, point_estimate)
    ubi <- c(point_estimate, point_estimate * 10)
  } else {
    lbi <- c(point_estimate * 10, point_estimate)
    ubi <- c(point_estimate, point_estimate / 10)
  }

  lb <- -Inf
  if (type == "interval" || type == "lower_bound") {
    lb <- stats::uniroot(
      f = cost,
      interval = lbi,
      extendInt = "upX"
    )$root
  }

  ub <- Inf
  if (type == "interval" || type == "upper_bound") {
    ub <- stats::uniroot(
      f = cost,
      interval = ubi,
      extendInt = "downX"
    )$root
  }

  c(lb, ub)
}

get_ci <- function(object,
                   pf,
                   conf_level = 0.95,
                   regression = FALSE,
                   ...) {
  if (!inherits(object, "param"))
    abort("The first argument should be a `param` object.")
  point_estimate <- object$point_estimate
  rngs <- dials::range_get(object, original = FALSE)
  if (!dials::is_unknown(rngs$lower) && !dials::is_unknown(rngs$upper))
    return(object)
  type <- if (!dials::is_unknown(rngs$lower))
    "upper_bound"
  else if (!dials::is_unknown(rngs$upper))
    "lower_bound"
  else
    "interval"
  ci <- compute_confidence_interval(
    pf = pf,
    conf_level = conf_level,
    type = type,
    point_estimate = point_estimate,
    regression = regression
  )
  if (!dials::is_unknown(rngs$lower))
    rngs$upper <- ci[2]
  else if (!dials::is_unknown(rngs$upper))
    rngs$lower <- ci[1]
  else {
    rngs$lower <- ci[1]
    rngs$upper <- ci[2]
  }
  res <- dials::range_set(object, rngs)
  res$point_estimate <- point_estimate
  class(res) <- c("inferred_param", class(res))
  res
}
