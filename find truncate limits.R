N <- 5000
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = N, clear = FALSE, width= 60)


# defining global parameters
n <- 100
beta0 <- 0
beta <- c(0.5, 1)
gamma0 <- 1
gamma <- c(1, 1)
phi <- 5
theta <- 0.4

#tuning parameters
rSq_theoretical <- 0.90  # used to generate sigma 
rSq_latent <- 0.60      # used to generate scale

sigma <- theory_rSq_to_sigma(beta, phi, theta, rSq_theoretical)

pr_bounds <- matrix(nrow = N, ncol = 2)
colnames(pr_bounds) <- c('lower', 'upper')

for (i in 1:N) {
  pb$tick()
  
  df <- generate_full_data(beta0, beta, phi, theta, sigma, n, misspecified = FALSE)
  df$z <- as.character(df$z)
  s <- latent_rSq_to_logistic_scale(gamma, rSq_latent, df)
  df.completeness <- generate_completeness(gamma0, gamma, df, s)
  df.ods <- generate_ods_data(df, df.completeness)
  pr_bound <- pr_observe(df.completeness, df.ods) %>%
    quantile(probs = c(0.025, 0.975))
  
  pr_bounds[i,1] <- pr_bound[1]
  pr_bounds[i,2] <- pr_bound[2]
}

pr_bounds <- as.data.frame(pr_bounds)

ggplot(data = pr_bounds) +
  geom_density(aes(x = lower), fill = I("red"), alpha = 0.5) +
  geom_density(aes(x = upper), fill = I("blue"), alpha = 0.5) +
  theme_bw()

pr_bounds %>%
  summarise(
    across(everything(), list(
      p025 = ~quantile(., 0.025),
      p975 = ~quantile(., 0.975))
      )
    )




