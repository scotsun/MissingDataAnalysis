library(tidyverse)

# TODO: update plot functions for the EM algorithm

mse <- function(data, beta_hat) {
  data <- data[complete.cases(data),]
  y_hat <- beta_hat[1] + beta_hat[2]*data$x + beta_hat[3]*as.numeric(data$z)
  return(mean((y_hat - data$y)^2))
}

plot_cc_vs_full <- function(df.ods, df) {
  cc <- lm(y ~ x + z, data = df.ods)
  cc.coef <- coefficients(cc)
  df <- df %>% mutate(z = paste('full', z))
  df.full_lm <- lm(y ~ x + z, df)
  p <- ggplot() +
    geom_point(data = na.omit(df.ods), aes(x = x, y = y, shape = z, color = z)) +
    geom_abline(intercept = cc.coef['(Intercept)'], 
                slope = cc.coef['x'],
                color="#F5B7B1", 
                linetype="solid", size=1.5) +
    geom_abline(intercept = cc.coef['(Intercept)'] + cc.coef['z1'], 
                slope = cc.coef['x'], 
                color="#AED6F1", 
                linetype="solid", size=1.5) +
    geom_smooth(data = df, 
                aes(x = x, y = predict(df.full_lm, df), group = z, color = z), linetype = 2,
                formula = 'y ~ x', method = "lm", se = FALSE, fullrange = TRUE) +
    scale_color_manual(values = c("#F5B7B1","#AED6F1", "red", "blue"), 
                       aesthetics = "color") +
    scale_shape(guide = 'none') +
    labs(color = '') +
    theme_bw() + 
    theme(legend.position = 'bottom')
  return(p)
}

plot_ipw_vs_full <- function(df.ods, df) {
  ipw.coef <- IPW(y ~ x + z, df.ods)
  df <- df %>% mutate(z = paste('full', z))
  df.full_lm <- lm(y ~ x + z, df)
  p <- ggplot() +
    geom_point(data = na.omit(df.ods), aes(x = x, y = y, shape = z, color = z)) +
    geom_abline(intercept = ipw.coef[1], 
                slope = ipw.coef[2],
                color="#F5B7B1", 
                linetype="solid", size=1.5) +
    geom_abline(intercept = ipw.coef[1] + ipw.coef[3], 
                slope = ipw.coef[2], 
                color="#AED6F1", 
                linetype="solid", size=1.5) +
    geom_smooth(data = df, 
                aes(x = x, y = predict(df.full_lm, df), group = z, color = z), linetype = 2,
                formula = 'y ~ x', method = "lm", se = FALSE, fullrange = TRUE) +
    scale_color_manual(values = c("#F5B7B1","#AED6F1", "red", "blue"), 
                       aesthetics = "color") +
    scale_shape(guide = 'none') +
    labs(color = '') +
    theme_bw() + 
    theme(legend.position = 'bottom')
  return(p)
}

plot_em_vs_full <- function(tolerance, df.ods, df, verbose) {
  em.coef <- EM(tolerance, df.ods, verbose)
  df <- df %>% mutate(z = paste('full', z))
  df.full_lm <- lm(y ~ x + z, df)
  p <- ggplot() +
    geom_point(data = na.omit(df.ods), aes(x = x, y = y, shape = z, color = z)) +
    geom_abline(intercept = em.coef[1], #TODO: consider to generalize in the future
                slope = em.coef[2],
                color="#F5B7B1", 
                linetype="solid", size=1.5) +
    geom_abline(intercept = em.coef[1] + em.coef[3], 
                slope = em.coef[2], 
                color="#AED6F1", 
                linetype="solid", size=1.5) +
    geom_smooth(data = df, 
                aes(x = x, y = predict(df.full_lm, df), group = z, color = z), linetype = 2,
                formula = 'y ~ x', method = "lm", se = FALSE, fullrange = TRUE) +
    scale_color_manual(values = c("#F5B7B1","#AED6F1", "red", "blue"), 
                       aesthetics = "color") +
    scale_shape(guide = 'none') +
    labs(color = '') +
    theme_bw() + 
    theme(legend.position = 'bottom')
  return(p)
}

get_plot_data <- function(parameter, simulation_df) {
  n <- dim(simulation_df)[1]
  simulation_df_subset <- simulation_df %>% 
    select(starts_with(parameter))
  plot_data <- data.frame(
    estimate = simulation_df_subset %>% 
      unlist() %>% 
      unname(),
    approach = c(rep("full", n), rep("cc", n), rep("ipw", n), rep("em", n)),
    iter = rep(seq(1,n), 4)
  )
  return(plot_data)
}

trace_plot_simulated_estimates <- function(parameter, simulation_df, alpha = 0.3,
                                           colorless_shadow = FALSE, smooth_method = "lm") {
  p <- ggplot(data = get_plot_data(parameter, simulation_df),
              aes(x = iter, y = estimate, color = approach))
  if (colorless_shadow) {
    p <- p + geom_line(alpha = alpha, color = "grey")
  } else {
    p <- p + geom_line(alpha = alpha)
  }
  if (smooth_method == "loess") {
    p <- p + geom_smooth(method = smooth_method, formula = 'y ~ x', se = FALSE, linetype = "dashed")
  } else if (smooth_method == "lm") {
    p <- p + geom_smooth(method = smooth_method, formula = 'y ~ 1', se = FALSE, linetype = "dashed")
  }
  p <- p + scale_colour_manual(values = c("red", "orange", "cornflowerblue", "chartreuse3")) +
    theme_bw()
  return(p)
}
