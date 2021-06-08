library(tidyverse)

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

plot_em_vs_full <- function(df.ods, df, ...) {
  em.coef <- EM(df.ods, ...)
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
    approach = c(rep("Full Data", n), rep("CC", n), rep("IPW", n), rep("ML (EM)", n)),
    iter = rep(seq(1,n), 4)
  )
  return(plot_data)
}


density_plot_simulated_estimates <- function(parameter, simulation_df, alpha = 0.5) {
  approach_colors <- c("CC" = "red", "ML (EM)" = "orange", 
                  "Full Data" = "cornflowerblue", "IPW" = "chartreuse3")
  approach_lty <- c("CC" = "dotted", "ML (EM)" = "longdash",
                    "Full Data" = "solid", "IPW" = "dashed")
  p <- ggplot(data = get_plot_data(parameter, simulation_df),
              aes(x = estimate)) +
    geom_density(aes(color = approach, fill = approach), alpha = alpha) +
    scale_color_manual(values = approach_colors) +
    scale_fill_manual(values = approach_colors)
  expected_est <- get_plot_data(parameter, simulation_df) %>% 
    group_by(approach) %>% 
    summarise(mean_est = mean(estimate))
  for (i in 1:nrow(expected_est)) {
    p <- p + geom_vline(xintercept = as.numeric(expected_est[i,2]), 
                        color = approach_colors[as.character(expected_est[i,1])],
                        linetype = approach_lty[as.character(expected_est[i,1])],
                        size = 1)
  }
  p <- p + theme_bw()
  return(p)
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
    p <- p + geom_smooth(aes(linetype = approach), method = smooth_method, formula = 'y ~ x', se = FALSE)
  } else if (smooth_method == "lm") {
    p <- p + geom_smooth(aes(linetype = approach), method = smooth_method, formula = 'y ~ 1', se = FALSE)
  }
  p <- p + scale_color_manual(values = c("red", "orange", "cornflowerblue", "chartreuse3")) +
    theme_bw()
  return(p)
}

coverage <- function(est_df, se_df, param) {
  N <- dim(est_df)[1]
  param_df <- data.frame(
    cbind(rep(param[1], N), rep(param[2], N), rep(param[3], N))
  )
  colnames(param_df) <- c("b0", "b1", "b2")
  l <- est_df - 1.96 * se_df
  u <- est_df + 1.96 * se_df
  ((param_df > l) & (param_df < u)) %>% apply(2, mean) 
}

MSE_estimator <- function(est_df, param) {
  expected_est <- est_df %>% apply(2, mean)
  bias <- expected_est - param
  var_est <- est_df %>% 
    var() %>% 
    unname() %>% 
    diag()
  return(bias^2 + var_est)
}
