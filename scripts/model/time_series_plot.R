
# add column for time series direction (left or right) 
# add starting and ending points for geom_segment i.e. to colour the different time series segments
ts_data$direction <- NA
ts_data$xend <- NA
for (i in 1:(nrow(ts_data)-1)) {
  if (ts_data$subject[i+1] == ts_data$subject[i]) {
    ts_data$xend[i] <- ts_data$state[i+1]
    
    if (ts_data$state[i+1] > ts_data$state[i]) {
      ts_data$direction[i] <- "right"
    } else {
      ts_data$direction[i] <- "left"
    }
  } else next
}

ts_data$mean_val <- NA
for (i in unique(ts_data$subject)) {
  subject_vals <- which(ts_data$subject == i)
  ts_data$mean_val[subject_vals] <- (min(ts_data$state[subject_vals]) + max(ts_data$state[subject_vals]))/2
}

ts_data$mean_val <- factor(ts_data$mean_val)

ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=subject, group=as.character(subject)), colour = "royalblue4") +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "State", y = "Subject") +
  theme_classic() + 
  #geom_vline(xintercept = real_root, linetype = "dotted") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# make plots
cols <- c("red4", "royalblue4", "turquoise2", "orchid4", "darkorange2", "yellow2")

ts_plot_ordered <- ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "royalblue4", size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# add tipping point line
#if (ground_truth) {
#  ts_plot_ordered <- ts_plot_ordered + geom_vline(xintercept = real_root, linetype = "dotted")
#}
