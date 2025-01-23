
# calculated based on drift and diffusion posterior draws
sd_draws <- data.frame("x" = x_pred)
# calculate stationary density for each drift and diffusion draw
#for (i in 1:(ncol(drift_draws)-1)) {
#  sd_draws[, i+1] <- nonparametric_stationary_density(drift_draws[, i], 
#                                                      diff_draws[, i], 
#                                                      x_pred)          
#}

# get stationary density
pos_density <- nonparametric_stationary_density_from_stan_samples_GP(x_pred, samples, pred_inf)
# or using drift and diffusion posterior means
stat_dens <- nonparametric_stationary_density(drift[,1], diffusion[,1], x_pred) 

# plot
density_plot <- ggplot(data = pos_density) +
  #geom_line(data = sd_draws[,2:ncol(sd_draws)] %>% 
  #            cbind(x = sd_draws$x) %>% 
  #            melt("x") %>%
  #            cbind(lp = rep(log_prob, each = length(x_pred))), 
  #          aes(x = x, y = value, group = variable, alpha = exp(lp)), size = 0.2, colour = "#f0d4ec") +
  geom_ribbon(aes(x = x_pred, ymin = lower_2.5, ymax = upper_97.5), fill = "thistle2", alpha = 0.7) + 
  geom_ribbon(aes(x = x_pred, ymin = lower_25, ymax = upper_75), fill = "#c6228f", alpha = 0.3) +
  labs(title = "", x = "System state", y = "Density") +
  geom_line(aes(x, mean), alpha = 1, size = 0.7, colour = "#c6228f", show.legend = FALSE) +
  #geom_line(aes(x_pred, stat_dens), alpha = 1, size = 0.7, colour = "#c6228f", show.legend = FALSE) +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 1.7)) +
  guides(alpha = FALSE) +
  theme_classic()

# plot effective potential
potential_plot <- ggplot() +
  labs(title = "", x = "State", y = "Potential") +
  #geom_line(data = sd_draws[,2:ncol(sd_draws)] %>% 
  #            cbind(x = sd_draws$x) %>% 
  #            melt("x"), 
  #          aes(x = x, y = -log(value), group = variable), size = 0.2, alpha = 0.4, colour = "#f0d4ec") +
  geom_line(data = pos_density, aes(x, -log(mean)), alpha = 1, size = 0.7, colour = "#c6228f", show.legend = FALSE) +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  theme_classic()

if (ground_truth) {
  density_plot <- density_plot + geom_line(data = ground_truth_stat_dens, aes(x, y), size = 0.7, linetype = "dashed")
  potential_plot <- potential_plot + geom_line(data = ground_truth_stat_dens, aes(x, -log(y)), size = 0.7, linetype = "dashed")
}
