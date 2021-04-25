get_inital_estimate <- function(data) {
  theta_hat <- unname((table(data$z) / sum(table(data$z)))['1'])
  sigma_hat <- sd(data$y)
  beta_hat <- unname(coef(lm(y ~ x + z, data))) # containing both beta0 and beta
  return(c(theta_hat, sigma_hat, beta_hat))
}

EM <- function(tolerance, data) {
  # augmented data used to update beta's
  augm_data <- augmented_data(data) 
  # since they are constant through out the process, so take it to the beginning
  n <- nrow(data)
  
  # preprocess the data and split them into 2 sets
  data$z <- as.numeric(data$z)
  data.R <- data[complete.cases(data),]
  data.R_bar <- data[!complete.cases(data),]
  num_cc <- nrow(data.R)
  
  estimate <- get_inital_estimate(data)
  phi_prime <- conditional_phi(estimate, data)
  Q_current <- 0
  Q_next <- Q(estimate, phi_prime, data.R, data.R_bar)
  iter <- 0
  while (abs(Q_current - Q_next) > tolerance) {
    iter <- iter + 1
    print(paste0("num of iter:", iter))
    Q_current <- Q_next
    
    A_ <- A(phi_prime, data.R)
    B_ <- B(phi_prime, data.R)
    S_ <- S(estimate, phi_prime, data.R, data.R_bar)
    
    n <- nrow(data)
    num_cc <- nrow(data.R)
    reg_coef <- lm(y ~ x + z, data = augm_data, weights = w(num_cc, phi_prime)) %>% 
      coef() %>% unname()
    
    estimate <- c(A_ / (A_ + B_), # theta
                  sqrt(S_/n),     # sigma
                  reg_coef)       # beta's
    
    print(estimate)
    
    phi_prime <- conditional_phi(estimate, data.R_bar)
    
    Q_next <- Q(estimate, phi_prime, data.R, data.R_bar)
  }
  return(estimate)
}

conditional_phi <- function(estimate, data) {
  data.R_bar <- data[!complete.cases(data), ]
  theta <- estimate[1]
  sigma <- estimate[2]
  beta0 <- estimate[3]
  beta1 <- estimate[4]
  beta2 <- estimate[5]
  f1 <- dnorm(data.R_bar$y, mean = beta0 + beta1 * data.R_bar$x + beta2, sd = sigma) * theta
  f2 <- dnorm(data.R_bar$y, mean = beta0 + beta1 * data.R_bar$x, sd = sigma) * (1 - theta)
  return(f1 / (f1 + f2))
}

Q <- function(estimate, phi_prime, data.R, data.R_bar) {
  theta <- estimate[1]
  sigma <- estimate[2]
  beta0 <- estimate[3]
  beta1 <- estimate[4]
  beta2 <- estimate[5]
  return(
    sum(dnorm(data.R$y, mean = beta0 + beta1 * data.R$x + beta2*as.integer(data.R$z), sd = sigma, log = TRUE)) +
      sum(dbinom(data.R$z, 1, prob = theta, log = TRUE)) +
      sum(dnorm(data.R_bar$y, mean = beta0 + beta1 * data.R_bar$x, sd = sigma, log = TRUE) * (1 - phi_prime)) +
      sum(dnorm(data.R_bar$y, mean = beta0 + beta1 * data.R_bar$x + beta2, sd = sigma, log = TRUE) * phi_prime) +
      sum((phi_prime * log(theta) + (1 - phi_prime) * log(1 - theta)))
  )
}

A <- function(phi_prime, data.R) {
  return(sum(data.R$z) + sum(phi_prime))
}

B <- function(phi_prime, data.R) {
  return(sum(1 - data.R$z) + sum(1 - phi_prime))
}

S <- function(estimate, phi_prime, data.R, data.R_bar) {
  beta0 <- estimate[3]
  beta1 <- estimate[4]
  beta2 <- estimate[5]
  return(
    sum((data.R$y - beta0 - beta1*data.R$x - beta2*data.R$z)^2) + 
      sum(((data.R_bar$y - beta0 - beta1*data.R_bar$x - beta2)*phi_prime)^2) +
      sum(((data.R_bar$y - beta0 - beta1*data.R_bar$x - beta2)*(1 - phi_prime))^2)
  )
}
