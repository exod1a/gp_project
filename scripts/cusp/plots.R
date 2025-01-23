
# create plots for Cusp model

# cusp time series plot
ts_plot <- ts_data %>%
  ggplot() +
  geom_line(aes(x = time, y = state), colour = "royalblue4") +
  #geom_point(aes(x = time, y = state), colour = "royalblue4") +
  labs(x = "Time", y = "System State") +
  geom_hline(yintercept = roots[2], linetype = "dotted", size = 0.5) +
  coord_cartesian(ylim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
  theme_classic()

# cusp drift plot
drift_plot <- ggplot() +
  geom_line(data = ground_truth_drift, aes(x, y), size = 0.7) + labs(title = "", x = "State", y = "Drift") + 
  geom_hline(yintercept=0, size = 0.4) +
  geom_point(aes(roots, 0), size = 2.5) +  
  geom_vline(xintercept = roots[2], linetype = "dotted", size = 0.5) +
  coord_cartesian(xlim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
  theme_classic()

# cusp diffusion plot
diff_plot <- ggplot() +
  geom_line(aes(ground_truth_diff$x, ground_truth_diff$y), size = 0.7) + labs(title = "", x = "State", y = "Diffusion") +
  coord_cartesian(xlim = c(min(ground_truth_diff$x), max(ground_truth_diff$x))) +
  geom_vline(xintercept = roots[2], linetype = "dotted", size = 0.5) +
  theme_classic()

# cusp density plot
dens_plot <- ggplot(data = ground_truth_stat_dens) +
  labs(title = "", x = "State", y = "Density") +
  geom_ribbon(aes(x, ymin = 0, ymax = y), fill = "#f0d4ec", alpha = 1, colour = "black") +
  coord_cartesian(xlim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
  geom_vline(xintercept = roots[2], linetype = "dotted", size = 0.5) +
  theme_classic()

# potential function plot
potential_plot <- ggplot() +
  labs(title = "", x = "State", y = "Potential") +
  geom_area(data = ground_truth_stat_dens, aes(x, -log(y)), alpha = 1, size = 0.5, 
            fill = "#f0d4ec", colour = "black", show.legend = FALSE) +
  geom_vline(xintercept = roots[2], linetype = "dotted", size = 0.5) +
  coord_cartesian(xlim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
  theme_classic()

# plot exit time
cusp_et_plot <- ggplot(data = ground_truth_et) +
  geom_line(aes(x, et), size = 0.7) +
  geom_segment(aes(x = ground_truth_et_r[[1]][1], y = mean_exit_time_right, 
                   xend = ground_truth_et_r[[1]][length(ground_truth_et_r[[1]])], 
                   yend = mean_exit_time_right), linetype = "dashed", colour = "black", size = 1) +
  geom_segment(aes(x = ground_truth_et_l[[1]][1], y = mean_exit_time_left, 
                   xend = ground_truth_et_l[[1]][length(ground_truth_et_l[[1]])], 
                   yend = mean_exit_time_left), linetype = "dashed", colour = "black", size = 1) +
  labs(title = "", x = "State", y = "Exit Time") +
  geom_vline(xintercept = roots[2], linetype = "dotted", size = 0.5) +
  coord_cartesian(xlim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
  theme_classic()

# create combined plot
cusp_plot <- ggarrange(ts_plot + coord_flip(ylim = c(min(ground_truth_drift$x), max(ground_truth_drift$x))) + 
                         theme(axis.text.x = element_blank(),
                               axis.title.x = element_blank()),
                       dens_plot + theme(axis.text.x = element_blank(),
                                         axis.title.x = element_blank()),
                       drift_plot + theme(axis.text.x = element_blank(),
                                          axis.title.x = element_blank()),
                       potential_plot + theme(axis.text.x = element_blank(),
                                              axis.title.x = element_blank()),
                       diff_plot, cusp_et_plot, ncol = 2, nrow = 3, align = "hv")
