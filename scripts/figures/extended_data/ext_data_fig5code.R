
# code to make extended data figure 5
# get upper and lower bounds on the exit time given all of the data to see how it compares with our prediction

tp_crossings <- function(ts, tp, index) {
  
  cross_to_left <- c()
  cross_to_right <- c()
  
  for (i in 1:(nrow(ts)-1)) {
    if (ts[i+1, index] < tp & ts[i, index] < tp) {
      next
    } else if (ts[i+1, index] > tp & ts[i, index] < tp) {
      cross_to_right <- c(cross_to_right, i+1)
    } else if (ts[i+1, index] > tp & ts[i, index] > tp) {
      next
    } else if (ts[i+1, index] < tp & ts[i, index] > tp) {
      cross_to_left <- c(cross_to_left, i+1)
    }
  }
  
  return(list(length(cross_to_left) + length(cross_to_right), 
              data.frame("crossings" = c(cross_to_left, cross_to_right), 
                         "direction" = c(rep("l", length(cross_to_left)), rep("r", length(cross_to_right))))))
}

chosen_taxa <- "Prevotella_melaninogenica_et_rel"
samples <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$samples
x_pred <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$x_pred

data_read_from_file <- 1
pred_inf <- 1
ground_truth <- 0
source(here("scripts", "model", "run_drift_diff.R"))

source(here("scripts", "model", "stability.R"))
if (length(pos_roots) > 0) {
  source(here("scripts", "model", "tipping_point.R"))
  source(here("scripts", "model", "exit_time.R"))
}

ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
# rename column
names(ts_data_full)[4] <- c("state")

# predicted tipping point location
tp_location <- mean_tp

# vectors for lower and upper bounds, including their x values
lower_bound <- c()
lb_x <- c()
upper_bound <- c()
ub_x <- c()

# isolate subjects
subjects <- unique(ts_data_full$subject)

# 4 for prevotella, 1 for cusp and mendota
ind <- 4

# go through each short time series 
for (i in subjects) {
  print(i)
  # subset one short time series at a time
  subsetted_df <- ts_data_full[which(ts_data_full$subject == i), ]
  
  # if it crosses the tipping point, look for upper bound
  if (tp_crossings(subsetted_df, tp_location, ind)[[1]] != 0) {
    # if it crosses only once can easily find the upper bound
    if (tp_crossings(subsetted_df, tp_location, ind)[[1]] == 1) {
      
      # index where it crosses
      crossing_point <- tp_crossings(subsetted_df, tp_location, ind)[[2]][1, 1]
      # calculate upper bound
      upper_bound <- c(upper_bound, subsetted_df$time[crossing_point])
      # find x value exit time was calculated at
      ub_x <- c(ub_x, subsetted_df$state[1])
    }
  } else {
    # if it doesn't cross the tipping point, look for lower bound
    lower_bound <- c(lower_bound, max(subsetted_df$time) - min(subsetted_df$time))
    # find x value exit time was calculated at
    lb_x <- c(lb_x, subsetted_df$state[1])
  }
} 

ts_data_full$mean_val <- NA
for (i in unique(ts_data_full$subject)) {
  subject_vals <- which(ts_data_full$subject == i)
  ts_data_full$mean_val[subject_vals] <- (min(ts_data_full$state[subject_vals]) + max(ts_data_full$state[subject_vals]))/2
}

ts_data_full$mean_val <- factor(ts_data_full$mean_val)


ts_plot_ordered <- ts_data_full %>%
  ggplot() +
  geom_line(aes(x = state, y = mean_val, group=as.character(mean_val)), colour = "royalblue4", size = 0.2) +
  geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.1) +
  labs(x = "System state", y = "Subject", title = chosen_taxa) +
  #geom_vline(xintercept = mean_tp, linetype = "dashed", colour = "orange") +
  theme_classic() + 
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


time_scales_plot <- ggplot() +
  geom_ribbon(aes(x = pos_exit_time_tp_CI$x[1:(length(x_pred) + 1)], ymin = ET_centred_pos[, 1], ymax = ET_centred_pos[, 3]), fill = "mediumpurple2", alpha = 0.5) + 
  geom_ribbon(aes(x = pos_exit_time_tp_CI$x[1:(length(x_pred) + 1)], ymin = ET_centred_pos[, 1], ymax = ET_centred_pos[, 2]), fill = "mediumpurple3", alpha = 0.5) + 
  geom_line(aes(c(et_l[[1]], et_r[[1]]), c(et_l[[2]], et_r[[2]])), colour = "mediumpurple", size = 0.7) +
  geom_point(aes(lb_x, lower_bound), colour = "snow3", alpha = 0.8) +
  geom_point(aes(ub_x, upper_bound), colour = "mediumturquoise", alpha = 0.8) +
  #geom_vline(xintercept = tp_location, linetype = "dashed", colour = "orange") +
  labs(title = "", x = "System state", y = "Exit time (days)") +
  coord_cartesian(ylim = c(0, 3000)) +
  theme_classic()

# create and save figure
figure <- ggarrange(ts_plot_ordered + coord_cartesian(xlim = c(min(x_pred), max(x_pred))),
                    time_scales_plot + coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 3000)),
                    ncol = 1, nrow = 2, align = "hv")

figure
# create path to save output
path <- here("output", "figures", "extended_data", "figure 5")
# save
ggsave(paste(path, "fig5.pdf", sep = "/"), figure, device = "pdf", height = 160, width = 110, units = "mm", dpi = 500)

