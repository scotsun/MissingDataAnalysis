library(tidyverse)
library(progress)

sigma_to_theory_rSq <- function(beta, phi, theta, sigma) {
  var_matrix <- matrix(c(phi^2, 0,
                             0, theta*(1-theta)), 
                       byrow = TRUE, nrow = 2, ncol = 2)
  var_explained <- t(beta) %*% var_matrix %*% beta
  return(var_explained / (var_explained + sigma^2))
}

theory_rSq_to_sigma <- function(beta, phi, theta, rSq) {
  var_matrix <- matrix(c(phi^2, 0,
                         0, theta*(1-theta)), 
                       byrow = TRUE, nrow = 2, ncol = 2)
  var_explained <- t(beta) %*% var_matrix %*% beta
  return(drop(sqrt(var_explained * (1 - rSq)/rSq)))
}

generate_full_data <- function(beta0, beta, phi, theta, sigma, n, misspecified = FALSE) {
  x <- rnorm(n, 0, phi)
  z <- rbinom(n, 1, theta)
  epsilon <- rnorm(n, 0, sigma)
  if (misspecified) {
    y <- beta0 + beta[1]*(sin(x) + x) + beta[2]*z + epsilon # mis-specified the mean function
  } else {
    y <- beta0 + beta[1]*x + beta[2]*z + epsilon
  }
  df <- data.frame(x = x, z = as.character(z), y = y)
  return(df)
}

logistic_scale_to_latent_rSq <- function(gamma, sigma, s, df) {
  var_matrix_lgr <- df %>% 
    select(x, y) %>% 
    cov()
  var_explained_lgr <- t(gamma) %*% var_matrix_lgr %*% gamma
  return(var_explained_lgr / (var_explained_lgr + 1/3*s^2*pi^2))
}

latent_rSq_to_logistic_scale <- function(gamma, rSq, df) {
  var_matrix_lgr <- df %>% 
    select(x, y) %>% cov() %>% unname()
  var_explained_lgr <- t(gamma) %*% var_matrix_lgr %*% gamma
  return(sqrt(3 * var_explained_lgr * (1 - rSq)/rSq)/pi)
}

generate_completeness <- function(gamma0, gamma, df, s) {
  n <- dim(df)[1]
  # latent variable (sigmoid + gamma0) in LGR
  sigmoid <- gamma[1]*df$x + gamma[2]*df$y + rlogis(n, location = 0, scale = s)
  completeness <- (sigmoid > -gamma0 + 0)
  return(completeness)
}

generate_ods_data <- function(df, completeness) {
  df <- df %>% mutate(z = ifelse(completeness == 1, z, NA))
  return(df)
}

###############################
# Inverse Probability Weighting
###############################

pr_observe <- function(completeness, data) {
  n <- dim(data)[1]
  lgr <- glm(completeness ~ data$y + data$x, family = "binomial", maxit=200)
  lgr.coef <- coef(lgr)
  data.mat <- c(rep(1, n), data$y, data$x) %>%
    matrix(nrow = nrow(data), ncol = ncol(data))
  return(exp(data.mat %*% lgr.coef) / (1 + exp(data.mat %*% lgr.coef)))
}

truncate <- function(pr, lower_bd, upper_bd) {
  for (i in 1:length(pr)) {
    if (pr[i] < lower_bd) {
      pr[i] <- lower_bd
    } else if (pr[i] > upper_bd) {
      pr[i] <- upper_bd
    } else {
      next
    }
  }
  return(pr)
}

IPW <- function(model, data, truncate_weight = FALSE, truncate_bounds = NULL) {
  pr <- pr_observe(as.integer(complete.cases(data)), data)
  if (truncate_weight) {
    pr <- truncate(pr, lower_bd = truncate_bounds[1], upper_bd = truncate_bounds[2])
  }
  ipw <- lm(y ~ x + z, data, 
            weights = 1/pr)
  ipw.coef <- ipw %>% 
    coefficients() %>% 
    unname()
  return(ipw.coef)
}

###############################
# IPW Inference
###############################

IPW_se <- function(model, data, truncate_weight = FALSE, truncate_bounds = NULL) {
  pr <- pr_observe(as.integer(complete.cases(data)), data)
  if (truncate_weight) {
    pr <- truncate(pr, lower_bd = truncate_bounds[1], upper_bd = truncate_bounds[2])
  }
  ipw <- lm(y ~ x + z, data, 
            weights = 1/pr)
  ipw.se <- ipw %>%
    vcov() %>%
    diag() %>%
    sqrt() %>% 
    unname()
  return(ipw.se)
}

###############################
# EM algorithm
###############################

get_initial_estimate <- function(data) {
  theta_hat <- unname((table(data$z) / sum(table(data$z)))['1'])
  sigma_hat <- sd(data$y)
  beta_hat <- unname(coef(lm(y ~ x + z, data))) # containing both beta0 and beta
  return(c(theta_hat, sigma_hat, beta_hat))
}

EM <- function(data, tolerance = 1e-6, max_iter = 1000, verbose = FALSE, report_allparam = FALSE) {
  # augmented data used to update beta's
  augm_data <- augmented_data(data) 
  # preprocess the data and split them into 2 sets
  data$z <- as.numeric(data$z)
  data.R <- data[complete.cases(data),]
  data.R_bar <- data[!complete.cases(data),]
  # since they are constant through out the process, so take it to the beginning
  n <- nrow(data)
  n_cc <- nrow(data.R)
  
  estimate <- get_initial_estimate(data)
  for (iter in 1:max_iter) {
    phi_prime <- conditional_phi(estimate, data)
    # to update theta
    A_ <- A(phi_prime, data.R)
    B_ <- B(phi_prime, data.R)
    # to update beta
    reg_coef <- lm(y ~ x + z, 
                   data = augm_data, 
                   weights = augment_weight(n_cc, phi_prime)) %>% 
      coef() %>% unname()
    # to update sigma
    S_ <- S(reg_coef, phi_prime, data.R, data.R_bar)
    #
    estimate_next <- c(A_ / (A_ + B_), # theta
                       sqrt(S_/n),     # sigma
                       reg_coef)       # beta's
    # check convergence
    check <- max(abs(estimate_next - estimate))
    if (check < tolerance) {
      estimate <- estimate_next
      break
    }
    phi_prime <- conditional_phi(estimate_next, data)
    estimate <- estimate_next
    # verbose print
    if (verbose) {
      print(paste0("num of iter:", iter))
      print(estimate)
    }
  }
  if (report_allparam) {
    return(list(estimate_theta = estimate[1], 
                estimate_sigma = estimate[2],
                estimate_beta = estimate[3:length(estimate)]))
  }
  return(estimate[3:length(estimate)])
}

conditional_phi <- function(estimate, data) {
  data.R_bar <- data[!complete.cases(data), ]
  
  theta <- estimate[1]
  sigma <- estimate[2]
  beta <- estimate[3:length(estimate)]
  
  mu0 <- drop(cbind(1, data.R_bar$x, 0) %*% beta)
  mu1 <- drop(cbind(1, data.R_bar$x, 1) %*% beta)
  
  f0 <- dnorm(data.R_bar$y, mean = mu0, sd = sigma)
  f1 <- dnorm(data.R_bar$y, mean = mu1, sd = sigma)
  return( (f1 * theta) / (f1 * theta + f0 * (1 - theta)))
}

A <- function(phi_prime, data.R) {
  return(sum(data.R$z) + sum(phi_prime))
}

B <- function(phi_prime, data.R) {
  return(sum(1 - data.R$z) + sum(1 - phi_prime))
}

S <- function(estimate_beta, phi_prime, data.R, data.R_bar) {
  X.R <- cbind(1, data.R$x, data.R$z)
  X.R_bar <- cbind(1, data.R_bar$x)
  return(
    sum((data.R$y - X.R %*% estimate_beta)^2) + 
      sum((data.R_bar$y - cbind(X.R_bar, 1) %*% estimate_beta)^2 * phi_prime) +
      sum((data.R_bar$y - cbind(X.R_bar, 0) %*% estimate_beta)^2 * (1 - phi_prime))
  )
}

augmented_data <- function(data) {
  data.R <- data[complete.cases(data),]
  data.R_bar_1 <- data[!complete.cases(data),]
  data.R_bar_0 <- data[!complete.cases(data),]
  
  data.R_bar_1$z <- '1'
  data.R_bar_0$z <- '0'
  
  return(rbind(data.R, data.R_bar_1, data.R_bar_0))
}

augment_weight <- function(num_cc, phi_prime) {
  return(
    c(rep(1, num_cc), phi_prime, 1 - phi_prime)
    )
}

###############################
# EM inference
###############################

f_j <- function(data, j, estimate_beta, estimate_sigma) {
  y <- data$y
  X <- cbind(1, data$x, j)
  y.mu <- drop(X %*% estimate_beta)
  return(dnorm(y, mean = y.mu, sd = estimate_sigma))
}

f_j_partial_sigma2 <- function(data, j, estimate_beta, estimate_sigma) {
  y <- data$y
  X <- cbind(1, data$x, j)
  y.mu <- drop(X %*% estimate_beta)
  f_j_ <- dnorm(y, mean = y.mu, sd = estimate_sigma)
  return(
    (-1/(2*estimate_sigma^2) + (y - y.mu)^2/(2*estimate_sigma^4)) * f_j_
  )
}

f_j_partial_beta <- function(data, j, estimate_beta, estimate_sigma) {
  y <- data$y
  X <- cbind(1, data$x, j)
  y.mu <- drop(X %*% estimate_beta)
  f_j_ <- dnorm(y, mean = y.mu, sd = estimate_sigma)
  return(
    (1/estimate_sigma^2 * (y - y.mu) * f_j_) * X
  )
}

fisher_score <- function(data, estimate) {
  estimate_theta <- estimate$estimate_theta
  estimate_sigma <- estimate$estimate_sigma
  estimate_beta <- estimate$estimate_beta
  z <- as.integer(data$z)
  is_complete <- complete.cases(data)
  
  f0 <- f_j(data, 0, estimate_beta, estimate_sigma)
  f1 <- f_j(data, 1, estimate_beta, estimate_sigma)
  
  D <- f1 - f0
  g <- f0 * (1 - estimate_theta) + f1 *  estimate_theta
  
  f0_partial_sigma2 <- f_j_partial_sigma2(data, 0, estimate_beta, estimate_sigma)
  f1_partial_sigma2 <- f_j_partial_sigma2(data, 1, estimate_beta, estimate_sigma)
  log_g_partial_sigma2 <- ((1 - estimate_theta) * f0_partial_sigma2 + estimate_theta * f1_partial_sigma2)/g
  
  f0_partial_beta <- f_j_partial_beta(data, 0, estimate_beta, estimate_sigma)
  f1_partial_beta <- f_j_partial_beta(data, 1, estimate_beta, estimate_sigma)
  log_g_partial_beta <- ((1 - estimate_theta) * f0_partial_beta + estimate_theta * f1_partial_beta)/g
  
  # initialization
  partial_theta <- D/g
  partial_sigma2 <- log_g_partial_sigma2
  partial_beta <- log_g_partial_beta
  # summation
  for (i in 1:length(is_complete)) {
    if (is_complete[i]) {
      partial_theta[i] <- partial_theta[i] + 
        z[i]*(1/estimate_theta - D[i]/g[i]) + 
        (1 - z[i])*(-1/(1 - estimate_theta) - D[i]/g[i])
      
      partial_sigma2[i] <- partial_sigma2[i] + 
        z[i]*(f1_partial_sigma2[i]/f1[i] - log_g_partial_sigma2[i]) + 
        (1 - z[i])*(f0_partial_sigma2[i]/f0[i] - log_g_partial_sigma2[i])
      
      partial_beta[i,] <- partial_beta[i,] + 
        z[i]*(f1_partial_beta[i,]/f1[i] - log_g_partial_beta[i,]) +
        (1 - z[i])*(f0_partial_beta[i,]/f0[i] - log_g_partial_beta[i,])
    } else {
      next
    }
  }
  score <- cbind(partial_theta, partial_sigma2, partial_beta)
  colnames(score) <- c("theta", "sigma2", "beta0", "beta1", "beta2")
  return(score)
}

EM_se <- function(data, estimate, report_allparam = FALSE) {
  n <- nrow(data)
  fisher_score_ <- fisher_score(data, estimate)
  obs_fisher_info <- n * var(fisher_score_)
  EM_vcov <- solve(obs_fisher_info)
  
  EM_all_se <- EM_vcov %>%
    diag() %>%
    sqrt()
  
  if (report_allparam) {
    return(EM_all_se)
  }
  
  idx <- grep('beta', names(EM_all_se), value = FALSE)
  return(unname(EM_all_se[idx]))
}

###############################
# Bootstrap
###############################

bootstrap <- function(B, model, data, method, boot_verbose, ...) {
  if (!method %in% c("cc", "ipw", "em")) {
    stop("undefined method")
  }
  n <- dim(data)[1]
  boot_estimates <- matrix(nrow = B, ncol = 3)
  if (boot_verbose) {
    pb <- progress_bar$new(
      format = "  downloading [:bar] :percent eta: :eta",
      total = B, clear = FALSE, width= 60)
  }
  for (b in 1:B) {
    if (boot_verbose) {pb$tick()}
    boot_data <- boot_sample(data)
    
    if (method == "cc") {
      estimates <- unname(coef(lm(model, boot_data)))
    } else if (method == "ipw") {
      # estimates <- IPW(model, boot_data)
      ################# TODO: debug
      tryCatch(
        expr = {estimates <- IPW(model, boot_data, ...)},
        error = function(e) {message("unknown source of error")}
      )
      #################
    } else if (method == "em") {
      estimates <- EM(boot_data, ...)
    } else {
      estimates <- rep(NA, 3)
    }
    boot_estimates[b,] <- estimates
  }
  return(boot_estimates)
}

boot_sample <- function(data) {
  data.complete <- data[complete.cases(data),]
  data.incomplete <- data[!complete.cases(data),]
  
  n.complete <- dim(data.complete)[1]
  n.incomplete <- dim(data.incomplete)[1]
  
  data.resample <- rbind(
    sample_n(data.complete, size = n.complete, replace = TRUE),
    sample_n(data.incomplete, size = n.incomplete, replace = TRUE)
  )
  
  data.resample <- data.resample %>%
    sample_n(size = n.complete + n.incomplete,
             replace = FALSE)
  if (nlevels(factor(data.resample$z)) < 2) {
    return(boot_sample(data))
  }
  return(data.resample)
}