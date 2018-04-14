#' @title CV for generalized supervised PCA
#'
#' @description
#' Run cross validation on dimension and \code{alpha} for generalized supervised PCA
#'
#' @param x covariate matrix
#' @param y response vector
#' @param ks the different dimensions \code{k} to try
#' @param alphas the different approximations to the saturated model \code{m} to try
#' @param family_x exponential family distribution of covariates
#' @param family_y exponential family distribution of response
#' @param folds if \code{folds} is a scalar, then it is the number of folds. If
#'  it is a vector, it should be the same length as the number of rows in \code{x}
#' @param quiet logical; whether the function should display progress
#' @param ... Additional arguments passed to \code{logisticPCA}
#'
#' @return A matrix of the CV deviance for the response with \code{k} in rows and
#'  \code{alpha} in columns
#'
#' @examples
#' # construct a low rank matrix in the logit scale
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_logit = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a binary matrix
#' mat = (matrix(runif(rows * cols), rows, cols) <= inv.logit.mat(mat_logit)) * 1.0
#'
#' \dontrun{
#' negloglikes = cv.gspca(mat, ks = 1:9)
#' plot(negloglikes)
#' }
#' @export
#' @importFrom stats predict
cv.gspca <- function(x, y, ks, alphas = 10^seq(-3, 5, by = .5),
                     family_x, family_y, folds = 5, quiet = TRUE, ...) {
  if (missing(family_x)) {
    family_x = exp_fam_guess(x)
  } else {
    check_family(x, family_x)
  }
  if (missing(family_x)) {
    family_y = exp_fam_guess(y)
  } else {
    check_family(y, family_y)
  }

  stopifnot(length(y) == nrow(x))

  if (length(folds) > 1) {
    # does this work if factor?
    if (length(unique(folds)) <= 1) {
      stop("If inputing CV split, must be more than one level")
    }
    if (length(folds) != nrow(x)) {
      stop("if folds is a vector, it should be of same length as nrow(x)")
    }
    cv = folds
  } else {
    stopifnot(folds >= 2 & folds <= (nrow(x) / 2))
    # TODO: stratify by non-missing response?
    cv = sample(1:folds, nrow(x), replace = TRUE)
  }

  deviances = matrix(0, length(ks), length(alphas),
                     dimnames = list(k = ks, alpha = alphas))
  for (k in ks) {
    for (alpha in alphas) {
      if (!quiet) {
        cat("k =", k, "alpha =", alpha, "")
      }
      for (c in unique(cv)) {
        if (!quiet) {
          cat(".")
        }
        gspca = genSupPCA(x[c != cv, ], y = y[c != cv], k = k, alpha = alpha,
                          family_x = family_x, family_y = family_y, ...)
        pred_theta = predict(gspca, newdat = x[c == cv, ], type = "link")
        deviances[k == ks, alpha == alphas] = deviances[k == ks, alpha == alphas] +
          exp_fam_deviance(x = y[c == cv], theta = pred_theta, family = family_y)
      }
      if (!quiet) {
        cat("", deviances[k == ks, alpha == alphas], "\n")
      }
    }
  }
  class(deviances) <- c("matrix", "cv.gspca")
  which_min = which(deviances == min(deviances), arr.ind = TRUE)
  if (!quiet) {
    cat("Best: k =", ks[which_min[1]], "alpha =", alphas[which_min[2]], "\n")
  }

  return(deviances)
}


#' @title Plot CV for generalized supervised PCA
#'
#' @description
#' Plot cross validation results generalized supervised PCA
#'
#' @param x a \code{cv.gspca} object
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrix in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#'
#' \dontrun{
#' loglikes = cv.gpca(mat, ks = 1:9, Ms = 3:6, family = "poisson")
#' plot(loglikes)
#' }
#' @export
#' @importFrom utils type.convert
plot.cv.gspca <- function(x, ...) {
  # replaces reshape2::melt(-x, value.name = "NegLogLikelihood")
  alphas = type.convert(colnames(x))
  ks = type.convert(rownames(x))
  df = data.frame(k = rep(ks, times = length(alphas)),
                  alpha = rep(alphas, each = length(ks)),
                  NegLogLikelihood = as.vector(x))

  if (ncol(x) == 1) {
    df$alpha = factor(df$alpha)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("k", "NegLogLikelihood", colour = "alpha")) +
      ggplot2::geom_line()
  } else {
    df$k = factor(df$k)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("alpha", "NegLogLikelihood", colour = "k")) +
      ggplot2::geom_line()
  }
  return(p)
}
