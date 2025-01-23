
# list of roots of the drift with a positive slope
tps <- list()

# take only roots for n_modes draws with exactly n_modes + 1 roots
num_roots <- num_roots[which(multistability == n_states)]

# agglomerate all positive roots
for (i in 1:length(num_roots)) {
  tps[[i]] <- x_pred[num_roots[[i]]$root[which(num_roots[[i]]$sign == 1)]]
}

# most probable tipping point
most_prob_tp <- mode_qi(unlist(tps))$y

# default setting .width = c(0.66, 0.95)
# make tipping point plot for the margin of the drift plot
CI_plot <- axis_canvas(drift_plot, axis = "x") +
           stat_pointinterval(aes(x = unlist(tps)), .width = c(0.50, 0.95), point_interval = mode_qi) +
           coord_cartesian(xlim = c(min(x_pred), max(x_pred)))

#####################################################################################################################
# find most probable tipping point
#tp_probabilities <- sort(table(unlist(tps))/length(unlist(tps))*100, decreasing = T)
#most_prob_tp <- as.numeric(names(tp_probabilities)[1])

#dens <- density(unlist(tps))
#df <- data.frame(x = dens$x, y = dens$y)

# supposedly most probable
#df$x[which(df$y == max(df$y))]
# median tp
#median_tp <- df$x[which(abs(df$x - median(unlist(tps))) == min(abs(df$x - median(unlist(tps)))))]
#mean_tp_2 <- mean(unlist(tps))

#line_data <- data.frame(xintercept = c(most_prob_tp, mean_tp, median_tp, mean_tp_2), 
#                        Lines = c("Most probable", "Drift posterior mean", "Median", "Mean"),
#                        color = c("blue", "red", "green", "purple"), stringsAsFactors = FALSE)


# with density plot
#dens_plot <- ggplot() +
#  geom_density(aes(unlist(tps))) +
#  geom_vline(data = line_data, aes(xintercept = xintercept, colour = Lines)) +
#  scale_colour_manual(values = line_data$color) +
 # coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
#  theme_classic()
  
#dens_plot <- insert_xaxis_grob(dens_plot, CI_plot, grid::unit(.06, "null"), position = "top")
#ggdraw(dens_plot)

#CI_plot <- ggplot() +
#           stat_pointinterval(aes(x = unlist(tps), y = rep(1, length(unlist(tps))))) +
#           coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(1, 10)) +
#           labs(title = "Tipping point region") +
#           theme_classic() +
#           easy_remove_axes()
