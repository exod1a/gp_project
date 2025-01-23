
# Fit model
if (!data_read_from_file) {
  # possible option for sampling statement "control = list(adapt_delta = 0.99)"
  # find how many cores available on system
  n_cores <- detectCores()
  # set number of cores to be the same as the number of chains or the maximum number of cores is below n_chains
  n_cores <- min(n_chains, n_cores)
  samples <- sampling(EM_model, stan_data, chains = n_chains, iter = n_iter, cores = n_cores, seed = 234)
}

# save
#saveRDS(samples, file = paste(path, "EM_data.rds", sep = "/"))

# get posterior mean
drift <- summary(samples, pars = "f_pred")[[1]][, c(1, 4:8)]

# get posterior draws
drift_draws <- rstan::extract(samples, "f_pred")[[1]]
drift_draws <- drift_draws %>% 
  t %>%
  as.data.frame() %>% 
  mutate(x = x_pred)

# find positive and negative roots of posterior mean
pos_roots <- drift_roots(drift[, 1], x_pred)$root[which(drift_roots(drift[, 1], x_pred)$sign == 1)]
neg_roots <- drift_roots(drift[, 1], x_pred)$root[which(drift_roots(drift[, 1], x_pred)$sign == -1)]

# tipping point of drift posterior mean
mean_tp <- x_pred[pos_roots]
log_prob <- rstan::extract(samples, "lp__")[[1]]

# plot drift
drift_plot <- ggplot() +
  #geom_line(data = drift_draws[,] %>% 
  #            melt("x") %>% 
  #            cbind(lp = rep(log_prob, each = length(x_pred))), 
  #          aes(x = x, y = value, group = variable, alpha = (exp(lp - (max(log_prob))) / (1 + exp(lp - (max(log_prob)))))), size = 0.2, colour = "#7a86c0") +
  geom_ribbon(aes(x = x_pred, ymin = drift[, 2], ymax = drift[, 6]), fill = "cyan", alpha = 0.5) +
  geom_ribbon(aes(x = x_pred, ymin = drift[, 3], ymax = drift[, 5]), fill = "blue", alpha = 0.5) + 
  geom_line(aes(x_pred, drift[, 1]), colour = "royalblue3", size = 1) +
  labs(x = "System state", y = "Drift") +
  geom_hline(yintercept = 0) +
  #geom_point(aes(ts_data$state, ts_data$drift), size = 2, alpha = 0.7, shape = 1) +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  guides(alpha = FALSE) +
  theme_classic()

# same for diffusion
# get posterior mean
diffusion <- 0.5*transform_diff(samples, pred_inf)
# get posterior draws
diff_draws <- 0.5*exp(rstan::extract(samples, "g_pred")[[1]])
diff_draws <- diff_draws %>% 
  t %>%
  as.data.frame() %>% 
  mutate(x = x_pred)

diff_plot <- ggplot() +
  labs(title = "", x = "System state", y = "Diffusion") +
  #geom_line(data = diff_draws[,] %>% 
  #            melt("x") %>% 
  #            cbind(lp = rep(log_prob, each = length(x_pred))), 
  #          aes(x = x, y = value, group = variable, alpha = (exp(lp - (max(log_prob))) / (1 + exp(lp - (max(log_prob)))))), size = 0.2, colour = "pink2") +
  geom_ribbon(data = diffusion, aes(x = x_pred, ymin = lower_2.5, ymax = upper_97.5), fill = "#90d9b0", alpha = 0.5) +
  geom_ribbon(data = diffusion, aes(x = x_pred, ymin = lower_25, ymax = upper_75), fill = "#4cc17c", alpha = 0.4) + 
  geom_line(data = diffusion, aes(x_pred, mean), colour = "#3dbc75", size = 0.7) + ##c6228f
  #geom_point(aes(ts_data$state, ts_data$diff), size = 2, alpha = 0.7, shape = 1) +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  guides(alpha = FALSE) +
  theme_classic()

# add ground truth to plots
if (ground_truth) {
  drift_plot <- drift_plot + geom_line(data = ground_truth_drift, aes(x, y), size = 0.7, linetype = "dashed")
  diff_plot <- diff_plot + geom_line(data = ground_truth_diff, aes(x, y), size = 0.7, linetype = "dashed")
}

if (!identical(pos_roots, numeric(0))) {
  drift_plot <- drift_plot + geom_point(aes(x_pred[pos_roots], drift[pos_roots, 1]), size = 2, colour = "black", fill = "white", shape = 21)
}
if (!identical(neg_roots, numeric(0))) {
  drift_plot <- drift_plot + geom_point(aes(x_pred[neg_roots], drift[neg_roots, 1]), size = 2, colour = "black", fill = "black", shape = 21)
}
