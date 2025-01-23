
# Need to solve the differential equation (Boundary Value Problem) D1(x)∂T/∂x + D2(x)∂^2T/∂x^2 = -1
# where D1(x) is the drift and D2(x) is the diffusion (we are now defining it the same way as them)

# for many posterior draws
# isolate posterior draws with highest probability number of modes
drift_mode_draws <- as.data.frame(cbind(drift_draws[, c(which(multistability == n_states))], x_pred))
#drift_mode_draws <- drift_mode_draws[, -to_remove]
diff_mode_draws  <- as.data.frame(cbind(diff_draws[, c(which(multistability == n_states))], x_pred))
#diff_mode_draws  <- diff_mode_draws[, -to_remove]

# run over all draws
# exit time posterior
pos_exit_time <- c("x" = NA, "ET" = NA, "draw" = NA)
pos_mean_exit_time <- list()
for (i in 1:(ncol(drift_mode_draws)-1)) {
  
  # find x values in each basin
  x_basins <- list()
  # find functions for differential equation
  p_basins <- list()
  r_basins <- list()
  exit_time <- list()

  # four cases for finding widths of basins:
  # 1 - first root is negative, last root is negative
  if (num_roots[[i]]$sign[1] == -1 && num_roots[[i]]$sign[nrow(num_roots[[i]])] == -1 && nrow(num_roots[[i]]) > 1) {
    
    # find locations of positive roots
    pos_equilibria <- num_roots[[i]]$root[which(num_roots[[i]]$sign == 1)]
    
    # need to include try functions because some posterior draws may lead to singular matrices
    # which are uninvertible when trying to calculate the exit time. This usually occurs once, twice,
    # or not at all for 4000 draws
    # left edge
    val <- try(exit_time_left_edge(x_pred, pos_equilibria, drift_mode_draws[, i], diff_mode_draws[, i]))
    if (inherits(val, "try-error")) {
      next
    }
    x_basins[[1]]  <- val[[1]]
    exit_time[[1]] <- val[[2]]
    
    # intermediate basins
    val <- try(exit_time_intermediate(x_pred, pos_equilibria, drift_mode_draws[, i], diff_mode_draws[, i]))
    if (inherits(val, "try-error")) {
      next
    }
    x_basins  <- list(x_basins, val[[1]])
    exit_time <- list(exit_time, val[[2]])
    
    # right edge
    val <- try(exit_time_right_edge(x_pred, pos_equilibria, drift_mode_draws[, i], diff_mode_draws[, i]))
    if (inherits(val, "try-error")) {
      next
    }
    x_basins[[length(x_basins) + 1]]   <- val[[1]]
    exit_time[[length(exit_time) + 1]] <- val[[2]]
    
    # save as data frame
    pos_exit_time <- rbind(pos_exit_time, 
                           data.frame("x" = unlist(x_basins), 
                                      "ET" = unlist(exit_time), 
                                      "draw" = rep(as.character(colnames(drift_mode_draws)[i]), length(unlist(exit_time)))))
  } else next
}
# remove first row
pos_exit_time <- pos_exit_time[-1, ]

# add log probability to exit time
# remove non-bistable draws
log_prob_bistable <- log_prob[which(multistability == n_states)]
# remove those that don't have exactly three roots
#log_prob_bistable <- log_prob_bistable[-to_remove]

# remove negative exit times since they are unphysical (if there are any)
if (sum(pos_exit_time$ET < 0)) {
  log_prob_bistable <- log_prob_bistable[-which(unique(pos_exit_time$draw) %in% unique(pos_exit_time$draw[which(pos_exit_time$ET < 0)]))]
  pos_exit_time <- pos_exit_time[-which((pos_exit_time$draw %in% unique(pos_exit_time$draw[which(pos_exit_time$ET < 0)]))), ]
}

# add to exit time data frame
pos_exit_time$lp <- rep(log_prob_bistable, each = length(x_pred) + 1)

# define tolerance range around mean tipping point within which to keep exit time posterior draws
tol = 0.1
to_keep <- which(unlist(tps) > mean_tp - tol & unlist(tps) < mean_tp + tol)
#unlist(tps)[to_keep]
drift_tp_CI_draws <- as.data.frame(cbind(drift_mode_draws[, to_keep], x_pred))
pos_exit_time_tp_CI <- pos_exit_time[which(pos_exit_time$draw %in% colnames(drift_tp_CI_draws)), ]

#percentage_to_take <- 0.01
#to_keep <- which(pos_exit_time$lp %in% sort(log_prob_bistable, decreasing = T)[1:floor(length(log_prob_bistable) * percentage_to_take)])
#pos_exit_time_reduced <- pos_exit_time[to_keep, ]

# calculate exit time for drift and diffusion posterior means
et_l <- exit_time_left_edge(x_pred, pos_roots, drift[,1], diffusion[,1])
et_r <- exit_time_right_edge(x_pred, pos_roots, drift[,1], diffusion[,1])

# plot exit time over all draws
et_pos_plot <- ggplot() +
               geom_line(data = pos_exit_time_tp_CI, aes(x = x, y = ET, group = draw), size = 0.2, colour = "steelblue2") +
               labs(title = "", x = "System state", y = "Exit time") +
               geom_line(aes(c(et_l[[1]], et_r[[1]]), c(et_l[[2]], et_r[[2]])), colour = "blue", size = 0.7) +
               coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 1000000)) +
               guides(alpha = FALSE) +
               theme_classic()
  
# exit time of cusp model
if (ground_truth) {
  et_pos_plot <- et_pos_plot + geom_line(data = ground_truth_et, aes(x, et), size = 0.7, linetype = "dashed")
}
et_pos_plot

upper_60 <- c(); upper_40 <- c(); lower_0 <- c()
leftover_draws <- length(unique(pos_exit_time_tp_CI$draw))

for (i in 1:(length(x_pred)+1)) {
  err           <- quantile(pos_exit_time_tp_CI[seq(from = i, to = leftover_draws*(length(x_pred)+1), by = length(x_pred)+1), 2], probs = c(0, 0.40, 0.60), na.rm = T)
  lower_0[i]    <- err[1]
  upper_40[i]   <- err[2]
  upper_60[i]   <- err[3]   
  
}
ET_centred_pos <- data.frame("lower_0" = lower_0, "upper_40" = upper_40, "upper_60" = upper_60)

# plot exit time over all draws
et_pos_plot <- ggplot() +
  geom_ribbon(aes(x = pos_exit_time_tp_CI$x[1:(length(x_pred) + 1)], ymin = ET_centred_pos[, 1], ymax = ET_centred_pos[, 3]), fill = "mediumpurple2", alpha = 0.5) + 
  geom_ribbon(aes(x = pos_exit_time_tp_CI$x[1:(length(x_pred) + 1)], ymin = ET_centred_pos[, 1], ymax = ET_centred_pos[, 2]), fill = "mediumpurple3", alpha = 0.5) + 
  labs(title = "", x = "System state", y = "Exit time") +
  geom_line(aes(c(et_l[[1]], et_r[[1]]), c(et_l[[2]], et_r[[2]])), colour = "mediumpurple", size = 0.7) +
  #geom_line(aes(x = pos_exit_time_tp_CI$x[1:(length(x_pred) + 1)], y = ET_centred_pos$mean), colour = "red", size = 0.7) +
  #coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 3000)) +
  guides(alpha = FALSE) +
  theme_classic()
et_pos_plot

#ggsave(paste(path, chosen_taxon, "exit_time.png", sep = "/"), et_pos_plot + geom_vline(xintercept = mean_tp, linetype = "dotted"), device = "png", height = 3, width = (7/4)*3, units = "in", dpi = 500)
# exit time of cusp model
if (ground_truth) {
  et_pos_plot <- et_pos_plot + geom_line(data = ground_truth_et, aes(x, et), size = 0.7, linetype = "dashed")
}
et_pos_plot
