#' Supervised dimensionality reduction for exponential family data
#'
#' @param x covariate matrix
#' @param y response vector
#' @param k dimension
#' @param alpha balance between dimensionality reduction of \code{x} and prediction of \code{y}
#' @param m value to approximate the saturated model in dimensionality reduction
#' @param family_x exponential family distribution of covariates
#' @param family_y exponential family distribution of response
#' @param quiet logical; whether the calculation should give feedback
#' @param max_iters maximum number of iterations
#' @param max_iters_per maximum iterations within each iteration
#' @param conv_criteria convergence criteria
#' @param discrete_deriv whether to calculate discrete derivatives w.r.t \code{U}
#'   instead of the closed form derivative with \code{beta} held constant
#' @param random_start whether to randomly initialize \code{U}
#' @param start_U initial value for \code{U}
#' @param start_beta initial value for \code{beta}
#' @param mu specific value for \code{mu}, the mean vector of \code{x}
#' @param grassmann logical; wether to optimize \code{U} on the Grassmann manifold.
#'   If \code{FALSE}, will optimize \code{U} on the Stiefel manifold
#'
#' @return An S3 object of class \code{gspca} which is a list with the
#' following components:
#' \item{mu}{the main effects for dimensionality reduction}
#' \item{U}{the \code{k}-dimentional orthonormal matrix with the loadings}
#' \item{beta}{the \code{k + 1} length vector of the coefficients}
#' \item{PCs}{the princial component scores}
#' \item{m}{the parameter inputed}
#' \item{family_x}{the exponential family of covariates}
#' \item{family_y}{the exponential family of response}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm.
#'    Should be non-increasing}
#'
#' @export
#' @importFrom RSpectra svds
#' @importFrom stats glm.fit gaussian binomial poisson
#' @import Matrix
#' @importFrom GrassmannOptim GrassmannOptim
#'
#' @examples
#' rows = 100
#' cols = 10
#' set.seed(1)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#'
#' # generate a count matrix and binary response
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#' response = rbinom(rows, 1, rowSums(mat) / max(rowSums(mat)))
#'
#' mod = genSupPCA(mat, response, k = 1, alpha = 0, family_x = "poisson", family_y = "binomial",
#'                 quiet = FALSE, max_iters_per = 5, grassmann = FALSE)
#'
#' plot(inv.logit.mat(cbind(1, mod$PCs) %*% mod$beta), response)
#' plot(rowSums(mat), response)
#' \dontrun{
#' ggplot(data.frame(PC = mod$PCs[, 1], y = response), aes(PC, y)) + stat_summary_bin(bins = 10)
#' }
genSupPCA <- function(x, y, k = 2, alpha = NULL, m = 4,
                           family_x = c("gaussian", "binomial", "poisson", "multinomial"),
                           family_y = c("gaussian", "binomial", "poisson"), quiet = TRUE,
                           max_iters = 100, max_iters_per = 5, conv_criteria = 1e-5, discrete_deriv = FALSE,
                           random_start = FALSE, start_U, mu, start_beta, grassmann = TRUE) {

  family_x = match.arg(family_x)
  family_y = match.arg(family_y)
  check_family(x, family_x)
  check_family(y, family_y)

  x = as.matrix(x)
  # missing_mat = is.na(x)
  n = nrow(x)
  d = ncol(x)
  ones = rep(1, n)

  stopifnot(NROW(y) == n)

  if (!grassmann && discrete_deriv) {
    warning("Cannot use discrete_deriv with Stiefel otpimization. Setting to FALSE")
    discrete_deriv = FALSE
  }

  # calculate the null log likelihood for % deviance explained and normalization
  null_theta_x = as.numeric(saturated_natural_parameters(matrix(colMeans(x, na.rm = TRUE), 1), family_x, m))
  null_theta_y = as.numeric(saturated_natural_parameters(matrix(mean(y, na.rm = TRUE), 1), family_y, Inf))

  # Initialize #
  ##################
  if (missing(mu)) {
    # eta = saturated_natural_parameters(x, family_x, m = m)
    # is.na(eta[is.na(x)]) <- TRUE
    # mu = colMeans(eta, na.rm = TRUE)
    mu = as.numeric(saturated_natural_parameters(matrix(colMeans(x, na.rm = TRUE), nrow = 1), family_x, Inf))
  }

  beta0 = saturated_natural_parameters(mean(y, na.rm = TRUE), family_y, Inf)

  if (is.null(alpha)) {
    null_loglike_x = exp_fam_deviance(x, outer(ones, mu), family_x)
    null_loglike_y = exp_fam_deviance(y, outer(ones, beta0), family_y)
    alpha = null_loglike_y / null_loglike_x
  }

  eta = saturated_natural_parameters(x, family_x, m = m) # + missing_mat * outer(ones, mu)
  eta_centered = scale(eta, center = mu, scale = FALSE)

  if (!missing(start_U)) {
    stopifnot(dim(start_U) == c(d, k))
    Qt = qr.Q(qr(start_U), complete = grassmann)
  } else if (random_start) {
    U = matrix(rnorm(d * k), d, k)
    Qt = qr.Q(qr(U), complete = grassmann)
  } else {
    if (grassmann) {
      udv = svd(eta_centered)
      Qt = udv$v
    } else {
      udv = RSpectra::svds(eta_centered, k)
      Qt = udv$v
    }
  }
  U = Qt[, 1:k, drop = FALSE]
  U_lag = U
  Qt_lag = matrix(0, d, d)

  if (missing(start_beta) | discrete_deriv) {
    beta = NULL
  } else {
    beta = start_beta
  }

  eta_sat_nat = saturated_natural_parameters(x, family_x, m = Inf)
  sat_loglike = exp_fam_log_like(x, eta_sat_nat, family_x)

  loss_trace = numeric(max_iters)
#   theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
#   loglike <- exp_fam_log_like(x, theta, family_x)
#   loss_trace[1] = (sat_loglike - loglike) / (n * d)
  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  if (grassmann) {
    params = list(
      Qt = Qt,
      dim = c(k, d),
      beta = beta,
      x = x,
      y = y,
      family_x = family_x,
      family_y = family_y,
      mu = mu,
      eta_centered = eta_centered,
      alpha = alpha
    )
  } else {
    tau = 1e-3
  }

  if (discrete_deriv) {
    max_iters_per = max(max_iters, max_iters_per)
    max_iters = 1
  }
  for (ii in seq_len(max_iters)) {
    if (!discrete_deriv) {
      mod_beta = suppressWarnings(
        stats::glm.fit(
          x = cbind(1, eta_centered %*% U),
          y = y,
          family = eval(parse(text = paste0("stats::", family_y, "()"))),
          control = list(maxit = max_iters_per),
          start = beta
        )
      )
      beta = mod_beta$coefficients
      if (grassmann) {
        params[["beta"]] <- beta
      }
    }

    if (grassmann) {
      mod_U = GrassmannOptim::GrassmannOptim(
        grassmann_objfun,
        params,
        max_iter = max_iters_per
      )

      if (is.null(mod_U$Qt)) {
        rnorm(1)
      }
      Qt = mod_U$Qt
      U = Qt[, 1:k, drop = FALSE]

      params[["Qt"]] <- Qt

      loss_trace[ii] <- -mod_U$fvalues[length(mod_U$fvalues)]
    } else {
      mod_U = FOptM::OptStiefelGBB(
        stiefel_objfun, U, opts = list(mxitr = max_iters_per, record = 0, tau = tau),
        x = x, y = y, family_x = family_x, family_y = family_y, mu = mu,
        eta_centered = eta_centered, beta = beta, alpha = alpha
      )

      U = mod_U$X
      tau = min(mod_U$out$tau * 10, 1e-3)

      loss_trace[ii] <- mod_U$out$fval
    }

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / ii * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(ii, "  ", loss_trace[ii], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }

    if (ii > 1 && sum((tcrossprod(U) %*% (diag(1, d, d) - tcrossprod(U_lag)))^2) < conv_criteria) {
      break
    } else {
      Qt_lag = Qt
      U_lag = U
    }
  }

  mod_beta = suppressWarnings(
    stats::glm.fit(
      x = cbind(1, eta_centered %*% U),
      y = y,
      family = eval(parse(text = paste0("stats::", family_y, "()"))),
      start = beta
    )
  )

  beta = mod_beta$coefficients

  if (grassmann) {
    mod_U$Qt <- NULL
    mod_U$fvalues <- -mod_U$fvalues
  } else {
    mod_U$out$fval <- mod_U$out$fval
  }

  object <- list(
    mu = mu,
    U = U,
    beta0 = beta0,
    beta = beta,
    PCs = eta_centered %*% U,
    theta_y = cbind(1, eta_centered %*% U) %*% beta,
    m = m,
    alpha = alpha,
    family_x = family_x,
    family_y = family_y,
    iters = ii,
    loss_trace = loss_trace[1:ii],
    Grassmann_object = mod_U
  )
  class(object) <- "gspca"
  object
}

#' @title Predict response with a generalized supervised PCA model
#'
#' @param object generalized supervised PCA object
#' @param newdata matrix of the same exponential family as covariates in \code{object}.
#'  If missing, will use the data that \code{object} was fit on
#' @param type the type of fitting required.
#'  \code{type = "link"} gives response variable on the natural parameter scale,
#'  \code{type = "response"} gives response variable on the response scale, and
#'  \code{type = "response"} gives matrix of principal components of \code{x}
#' @param ... Additional arguments
#' @examples
#' # construct a low rank matrices in the natural parameter space
#' rows = 100
#' cols = 10
#' set.seed(1)
#' loadings = rnorm(cols)
#' mat_np = outer(rnorm(rows), rnorm(cols))
#' mat_np_new = outer(rnorm(rows), loadings)
#'
#' # generate a count matrices and binary responses
#' mat = matrix(rpois(rows * cols, c(exp(mat_np))), rows, cols)
#' mat_new = matrix(rpois(rows * cols, c(exp(mat_np_new))), rows, cols)
#'
#' response = rbinom(rows, 1, rowSums(mat) / max(rowSums(mat)))
#' response_new = rbinom(rows, 1, rowSums(mat_new) / max(rowSums(mat_new)))
#'
#' mod = genSupPCA(mat, response, k = 1, family_x = "poisson", family_y = "binomial",
#'                 quiet = FALSE, max_iters_per = 1, discrete_deriv = FALSE)
#'
#' plot(predict(mod, mat, type = "response"), response_new)
#'
#' @export
predict.gspca <- function(object, newdata, type = c("link", "response", "PCs"), ...) {
  type = match.arg(type)

  if (missing(newdata)) {
    PCs = object$PCs
  } else {
    stopifnot(ncol(newdata) == nrow(object$U))
    check_family(newdata, object$family_x)

    eta = saturated_natural_parameters(newdata, object$family_x, m = object$m) # + missing_mat * outer(ones, mu)
    PCs = scale(eta, center = object$mu, scale = FALSE) %*% object$U
  }

  if (type == "PCs") {
    return(PCs)
  } else {
    theta = as.numeric(cbind(1, PCs) %*% object$beta)

    if (type == "link") {
      return(theta)
    } else if (type == "response") {
      return(exp_fam_mean(theta, object$family_y))
    }
  }
}
