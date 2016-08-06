#' Supervised matrix factorization for exponential family data
#'
#' @param x covariate matrix
#' @param y response vector
#' @param k dimension
#' @param alpha balance between dimensionality reduction of \code{x} and prediction of \code{y}
#' @param family_x exponential family distribution of covariates
#' @param family_y exponential family distribution of response
#' @param quiet logical; whether the calculation should give feedback
#' @param max_iters maximum number of iterations
#' @param conv_criteria convergence criteria
#' @param random_start whether to randomly initialize \code{A} and \code{B}
#' @param start_A initial value for \code{A}
#' @param start_B initial value for \code{B}
#' @param start_beta initial value for \code{beta}
#' @param mu specific value for \code{mu}, the mean vector of \code{x}
#'
#' @return An S3 object of class \code{gsmf} which is a list with the
#' following components:
#' \item{mu}{the main effects for dimensionality reduction}
#' \item{A}{the \code{n}x\code{k}-dimentional matrix with the scores}
#' \item{B}{the \code{d}x\code{k}-dimentional matrix with the loadings}
#' \item{beta}{the \code{k + 1} length vector of the coefficients}
#' \item{family_x}{the exponential family of covariates}
#' \item{family_y}{the exponential family of response}
#' \item{iters}{number of iterations required for convergence}
#' \item{loss_trace}{the trace of the average negative log likelihood of the algorithm.
#'    Should be non-increasing}
#'
#' @export
#' @importFrom RSpectra svds
#' @importFrom stats glm.fit gaussian binomial poisson
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
#' mod = genSupMF(mat, response, k = 1, alpha = 1, family_x = "poisson", family_y = "binomial",
#'                quiet = FALSE)
#'
#' plot(inv.logit.mat(cbind(1, mod$A) %*% mod$beta), response)
#' plot(rowSums(mat), response)
#' \dontrun{
#' ggplot(data.frame(PC = mod$PCs[, 1], y = response), aes(PC, y)) + stat_summary_bin(bins = 10)
#' }
genSupMF <- function(x, y, k = 2, alpha = NULL,
                      family_x = c("gaussian", "binomial", "poisson"),
                      family_y = c("gaussian", "binomial", "poisson"), quiet = TRUE,
                      max_iters = 1000, conv_criteria = 1e-5,
                      random_start = FALSE, start_A,  start_B, mu) {

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

  # Initialize #
  ##################
  if (missing(mu)) {
    mu = saturated_natural_parameters(colMeans(x, na.rm = TRUE), family_x, Inf)
  }

  beta0 = saturated_natural_parameters(mean(y, na.rm = TRUE), family_y, Inf)

  if (is.null(alpha)) {
    null_loglike_x = exp_fam_deviance(x, outer(ones, mu), family_x)
    null_loglike_y = exp_fam_deviance(y, outer(ones, beta0), family_y)
    alpha = null_loglike_y / null_loglike_x
  }

  if (!missing(start_B)) {
    stopifnot(dim(start_B) == c(d, k))
    B = start_B
  } else if (random_start) {
    B = matrix(rnorm(d * k), d, k)
  } else {
    udv = svd(scale(x, TRUE, FALSE))
    B = udv$v[, 1:k, drop = FALSE] %*% diag(sqrt(udv$d), k, k)
  }

  if (!missing(start_A)) {
    stopifnot(dim(start_A) == c(n, k))
    A = start_A
  } else if (random_start) {
    A = matrix(rnorm(n * k), n, k)
  } else {
    if (!missing(start_B)) {
      udv = svd(scale(x, TRUE, FALSE))
    }
    A = udv$u[, 1:k, drop = FALSE] %*% diag(sqrt(udv$d), k, k)
  }

  beta = NULL # as.numeric(solve(crossprod(cbind(1, A)), crossprod(cbind(1, A), y)))

  loss_trace = numeric(max_iters)
  #   theta = outer(ones, mu) + eta_centered %*% tcrossprod(U)
  #   loglike <- exp_fam_log_like(x, theta, family_x)
  #   loss_trace[1] = (sat_loglike - loglike) / (n * d)
  ptm <- proc.time()

  if (!quiet) {
    cat(0, "  ", loss_trace[1], "")
    cat("0 hours elapsed\n")
  }

  theta_x = outer(ones, mu) + tcrossprod(A, B)
  perfect_fit = FALSE

  for (ii in seq_len(max_iters)) {
    # update beta
    if (ii == 1) {
      mod_beta = suppressWarnings(
        stats::glm.fit(
          x = cbind(1, A),
          y = y,
          family = eval(parse(text = paste0("stats::", family_y, "()"))),
          control = list(maxit = 1),
          start = beta
        )
      )
      beta = mod_beta$coefficients
    } else {
      theta_y = cbind(1, A) %*% beta
      first_dir = exp_fam_mean(theta_y, family_y)
      second_dir = as.numeric(exp_fam_variance(theta_y, family_y))

      W = max(second_dir)
      Z = as.matrix(theta_y + (y - first_dir) / W)

      beta = as.numeric(solve(crossprod(cbind(1, A)) + diag(0.01, k + 1, k + 1),
                              crossprod(cbind(1, A), Z)))
    }

    # update A
    theta_y = cbind(1, A) %*% beta
    theta_x = outer(ones, mu) + tcrossprod(A, B)

    if (family_y == "binomial" && exp_fam_deviance(y, theta_y, family_y) < 1e-5) {
      loss_trace[ii] <- -sup_dim_red_log_like(x, y, theta_x, theta_y, family_x, family_y, alpha)
      break
    }

    first_dir_x = exp_fam_mean(theta_x, family_x)
    first_dir_y = exp_fam_mean(theta_y, family_y)
    second_dir_x = exp_fam_variance(theta_x, family_x)
    second_dir_y = exp_fam_variance(theta_y, family_y)

    if (alpha > 0) {
      theta = cbind(theta_x, theta_y)
      first_dir = cbind(first_dir_x, first_dir_y)
      second_dir = cbind(second_dir_x, second_dir_y)

      W = apply(second_dir, 2, max)
      W[-length(W)] <- W[-length(W)] * alpha
      Z = as.matrix(theta + (cbind(x, y) - first_dir) / outer(ones, W)) - outer(ones, c(mu, beta[1]))

      A = t(solve(t(rbind(B, beta[-1])) %*% diag(W) %*% rbind(B, beta[-1]) + diag(0.01, k, k),
                  t(rbind(B, beta[-1])) %*% diag(W) %*% t(Z)))
    } else {
      theta = theta_y
      first_dir = first_dir_y
      second_dir = second_dir_y

      W = max(second_dir)
      Z = as.matrix(theta + (y - first_dir) / W) - beta[1]

      A = t(solve(t(beta[-1]) %*% beta[-1] + diag(0.01, k, k),
                  t(beta[-1]) %*% t(Z)))
    }


    # update B
    first_dir = exp_fam_mean(theta_x, family_x)
    second_dir = exp_fam_variance(theta_x, family_x)

    W = apply(second_dir, 1, max)
    Z = as.matrix(theta_x + (x - first_dir) / outer(W, rep(1, d))) - outer(ones, mu)

    B = t(solve(t(A) %*% diag(W, n, n) %*% A + diag(0.01, k, k),
                t(A) %*% diag(W, n, n) %*% Z))

    # update A again
    theta_y = cbind(1, A) %*% beta
    theta_x = outer(ones, mu) + tcrossprod(A, B)

    first_dir_x = exp_fam_mean(theta_x, family_x)
    first_dir_y = exp_fam_mean(theta_y, family_y)
    second_dir_x = exp_fam_variance(theta_x, family_x)
    second_dir_y = exp_fam_variance(theta_y, family_y)

    if (family_y == "binomial" && -exp_fam_log_like(y, theta_y, family_y) < 1e-5) {
      perfect_fit = TRUE
    } else if (alpha > 0) {
      theta = cbind(theta_x, theta_y)
      first_dir = cbind(first_dir_x, first_dir_y)
      second_dir = cbind(second_dir_x, second_dir_y)

      W = apply(second_dir, 2, max)
      W[-length(W)] <- W[-length(W)] * alpha
      Z = as.matrix(theta + (cbind(x, y) - first_dir) / outer(ones, W)) - outer(ones, c(mu, beta[1]))

      A = t(solve(t(rbind(B, beta[-1])) %*% diag(W) %*% rbind(B, beta[-1]) + diag(0.01, k, k),
                  t(rbind(B, beta[-1])) %*% diag(W) %*% t(Z)))
    } else {
      theta = theta_y
      first_dir = first_dir_y
      second_dir = second_dir_y

      W = max(second_dir)
      Z = as.matrix(theta + (y - first_dir) / W) - beta[1]

      A = t(solve(t(beta[-1]) %*% beta[-1] + diag(0.01, k, k),
                  t(beta[-1]) %*% t(Z)))
    }

    # Calc Deviance
    theta_x = outer(ones, mu) + tcrossprod(A, B)
    theta_y = cbind(1, A) %*% beta

    loss_trace[ii] <- -sup_dim_red_log_like(x, y, theta_x, theta_y, family_x, family_y, alpha)

    if (!quiet) {
      time_elapsed = as.numeric(proc.time() - ptm)[3]
      tot_time = max_iters / ii * time_elapsed
      time_remain = tot_time - time_elapsed
      cat(ii, "  ", loss_trace[ii], "")
      cat(round(time_elapsed / 3600, 1), "hours elapsed. Max", round(time_remain / 3600, 1), "hours remain.\n")
    }

    if ((ii > 1 | perfect_fit) && abs(loss_trace[ii] - loss_trace[ii - 1]) < conv_criteria) {
      break
    } else {
      B_lag = B
    }
  }

  mod_beta = suppressWarnings(
    stats::glm.fit(
      x = cbind(1, A),
      y = y,
      family = eval(parse(text = paste0("stats::", family_y, "()"))),
      start = beta
    )
  )

  beta = mod_beta$coefficients

  object <- list(
    mu = mu,
    A = A,
    B = B,
    beta0 = beta0,
    beta = beta,
    theta_y = cbind(1, A) %*% beta,
    alpha = alpha,
    family_x = family_x,
    family_y = family_y,
    iters = ii,
    loss_trace = loss_trace[1:ii]
  )
  class(object) <- "gsmf"
  object
}

#' @title Predict response with a generalized supervised MF model
#'
#' @param object generalized supervised MF object
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
#' mat = matrix(rbinom(rows * cols, 1, c(inv.logit.mat(mat_np))), rows, cols)
#' mat_new = matrix(rbinom(rows * cols, 1, c(inv.logit.mat(mat_np_new))), rows, cols)
#'
#' response = rbinom(rows, 1, rowSums(mat) / max(rowSums(mat)))
#' response_new = rbinom(rows, 1, rowSums(mat_new) / max(rowSums(mat_new)))
#'
#' mod = genSupMF(mat, response, k = 2, alpha = 1000, family_x = "poisson", family_y = "binomial", quiet = FALSE)
#'
#' plot(predict(mod, type = "response"), response)
#' plot(predict(mod, mat_new, type = "response"), response_new)
#'
#' @export
predict.gsmf <- function(object, newdata, type = c("link", "response", "PCs"), quiet = TRUE,
                         max_iters = 1000, conv_criteria = 1e-5, start_A,...) {
  type = match.arg(type)

  if (missing(newdata)) {
    A = object$A
  } else {
    n = nrow(newdata)
    d = ncol(newdata)
    k = ncol(object$B)
    stopifnot(d == nrow(object$B))
    check_family(newdata, object$family_x)

    ones = rep(1, n)

    # solve for A
    last_deviance = Inf
    if (missing(start_A)) {
      theta = outer(ones, object$mu)
    } else {
      stopifnot(dim(start_A) == c(n, k))
      theta = outer(ones, object$mu) + tcrossprod(start_A, object$B)
    }

    for (ii in seq_len(max_iters)) {
      first_dir = exp_fam_mean(theta, object$family_x)
      second_dir = exp_fam_variance(theta, object$family_x)

      W = apply(second_dir, 2, max)
      Z = as.matrix(theta + (newdata - first_dir) / outer(ones, W)) - outer(ones, object$mu)

      A = t(solve(t(object$B) %*% diag(W, d, d) %*% object$B, t(object$B) %*% diag(W, d, d) %*% t(Z)))

      theta = outer(ones, object$mu) + tcrossprod(A, object$B)
      this_deviance = exp_fam_deviance(newdata, theta, object$family_x) / (n * d)

      if (!quiet) {
        cat(ii ," ", this_deviance, "\n")
      }

      if (abs(last_deviance - this_deviance) < 2 * conv_criteria)
        break
    }
  }

  if (type == "PCs") {
    return(A)
  } else {
    theta = as.numeric(cbind(1, A) %*% object$beta)

    if (type == "link") {
      return(theta)
    } else if (type == "response") {
      return(exp_fam_mean(theta, object$family_y))
    }
  }
}
