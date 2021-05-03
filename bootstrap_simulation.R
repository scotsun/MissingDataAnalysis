library(tidyverse)
library(progress)
source("./R/dataGenerator_ipw_em.R")

# defining global parameters
n <- 100
beta0 <- 0
beta <- c(0.5, 1)
gamma0 <- 1
gamma <- c(1, 1)
phi <- 5
theta <- 0.4

num_boots <- 800

#tuning parameters
rSq_theoretical <- 0.30  # used to generate sigma 
rSq_latent <- 0.6        # used to generate scale

N <- 500 #TODO: number of replication
seeds <- unlist(saved_seeds)

prefix <- c("b0_", "b1_", "b2_")
suffix <- c("se", "95lb", "95ub")
colnames_ <- c()
for (suf in suffix) {
  colnames_ <- append(colnames_, paste0(prefix, suf))
}


cc_bootstrap_simu <- matrix(nrow = N, ncol = length(colnames_))
ipw_bootstrap_simu <- matrix(nrow = N, ncol = length(colnames_))
em_bootstrap_simu <- matrix(nrow = N, ncol = length(colnames_))
colnames(cc_bootstrap_simu) <- colnames_
colnames(ipw_bootstrap_simu) <- colnames_
colnames(em_bootstrap_simu) <- colnames_

for (i in 1:N) {
  message(paste0(i, "/", N, " replications"))
  message("===========================================================")
  set.seed(seeds[i])
  sigma <- theory_rSq_to_sigma(beta, phi, theta, rSq_theoretical)
  df <- generate_full_data(beta0, beta, phi, theta, sigma, n, misspecified = FALSE)
  df$z <- as.character(df$z)
  s <- latent_rSq_to_logistic_scale(gamma, rSq_latent, df)
  df.completeness <- generate_completeness(gamma0, gamma, df, s)
  df.ods <- generate_ods_data(df, df.completeness)
  
  cc_boot_estimates <- bootstrap(num_boots, model = y ~ x + z, data = df.ods, 
                                 method = "cc", boot_verbose = TRUE)
  cc_bootstrap_simu[i,] <- cc_boot_estimates %>%
    as.data.frame() %>%
    summarise_all(.funs = list(se = ~ sd(x = .),
                               lb = ~ quantile(x = ., probs = 0.025),
                               ub = ~ quantile(x = ., probs = 0.975))) %>%
    unlist() %>% 
    unname()
  
  ipw_boot_estimates <- bootstrap(num_boots, model = y ~ x + z, data = df.ods, 
                                  method = "ipw", boot_verbose = TRUE,
                                  truncate_weight = TRUE,
                                  truncate_bounds = c(0.05, 0.99))
  ipw_bootstrap_simu[i,] <- ipw_boot_estimates %>%
    as.data.frame() %>%
    summarise_all(.funs = list(se = ~ sd(x = .),
                               lb = ~ quantile(x = ., probs = 0.025),
                               ub = ~ quantile(x = ., probs = 0.975))) %>%
    unlist() %>%
    unname()
  
  em_boot_estimates <- bootstrap(num_boots, model = y ~ x + z, data = df.ods, 
                                 method = "em", boot_verbose = TRUE,
                                 max_iter = 200)
  em_bootstrap_simu[i,] <- em_boot_estimates %>%
    as.data.frame() %>%
    summarise_all(.funs = list(se = ~ sd(x = .),
                               lb = ~ quantile(x = ., probs = 0.025),
                               ub = ~ quantile(x = ., probs = 0.975))) %>%
    unlist() %>%
    unname()
}

simulation_result <- list(seeds, cc_bootstrap_simu, em_bootstrap_simu, ipw_bootstrap_simu)


save(simulation_result, file = "./R/boot_simulation30.RData")


