# Grassman optimization ---------------------------------------------------

sup_dim_red_log_like <- function(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, U) {
  if (is.null(theta_y)) {
    mod_beta = suppressWarnings(
      stats::glm.fit(
        x = cbind(1, eta_centered %*% U),
        y = y,
        family = eval(parse(text = paste0("stats::", family_y, "()")))
      )
    )

    beta = mod_beta$coefficients
    theta_y = cbind(1, eta_centered %*% U) %*% beta
  }
  (-1 * exp_fam_deviance(y, theta_y, family_y) - alpha * exp_fam_deviance(x, theta_x, family_x)) #/ nrow(x) / (ncol(x) * alpha + 1)
}

sup_dim_red_directional_deriv <- function(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, beta, U, Z) {
  Ex = exp_fam_mean(theta_x, family_x)
  Ey = exp_fam_mean(theta_y, family_y)

  term1 = beta[-1] %*% crossprod(y - Ey, eta_centered) %*% Z
  term2a = crossprod(x - Ex, eta_centered)
  term2 = t(U) %*% (term2a + t(term2a)) %*% Z
  2 * (term1 + alpha * term2) / nrow(x) / (ncol(x) * alpha + 1)
}





grassmann_objfun <- function(W) {
  size = W$dim
  # print(size)
  k = size[1]
  d = size[2]
  if (is.null(W$Qt)) {
    rnorm(1)
  }
  U = W$Qt[, 1:k, drop = FALSE]
  Z = W$Qt[, (k + 1):d, drop = FALSE]
  theta_x = outer(rep(1, nrow(W$x)), W$mu) + W$eta_centered %*% tcrossprod(U)
  if (!is.null(W$beta)) {
    theta_y = cbind(1, W$eta_centered %*% U) %*% W$beta
  } else {
    theta_y = NULL
  }

  value <- sup_dim_red_log_like(W$x, W$y, theta_x, theta_y, W$family_x, W$family_y, W$alpha, W$eta_centered, U)
  if (!is.null(W$beta)) {
    gradient <- sup_dim_red_directional_deriv(W$x, W$y, theta_x, theta_y, W$family_x, W$family_y, W$alpha, W$eta_centered, W$beta, U, Z)
  }

  if (!is.null(W$beta)) {
    list(value = value,
         gradient = gradient)
  } else {
    list(value = value)
  }
}


# Stiefel optimization ----------------------------------------------------

sup_dim_red_deriv <- function(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, beta, U) {
  Ex = exp_fam_mean(theta_x, family_x)
  Ey = exp_fam_mean(theta_y, family_y)

  term1 = beta[-1] %*% crossprod(y - Ey, eta_centered)
  term2a = crossprod(x - Ex, eta_centered)
  term2 = t(U) %*% (term2a + t(term2a))
  -2 * t(term1 + alpha * term2) # / nrow(x) / (ncol(x) * alpha + 1)
}

stiefel_objfun <- function(U, x, y, family_x, family_y, mu, eta_centered, beta, alpha) {
  theta_x = outer(rep(1, nrow(x)), mu) + eta_centered %*% tcrossprod(U)
  theta_y = cbind(1, eta_centered %*% U) %*% beta

  value <- -sup_dim_red_log_like(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, U)
  gradient <- sup_dim_red_deriv(x, y, theta_x, theta_y, family_x, family_y, alpha, eta_centered, beta, U)

  list(f = value,
       G = gradient)
}
