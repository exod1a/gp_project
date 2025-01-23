
# Simulate from the cusp model
ts_data <- cusp_euler_maruyama(y0 = NULL,
                               times = time,
                               seed = NULL,
                               r = r,
                               alpha = alpha,
                               beta = beta,
                               lambda = lambda,
                               epsilon = epsilon) %>%
           cbind(time = time) %>%
           as.data.frame() %>% 
           set_colnames(c("state", "time"))

# subset data (otherwise it will take too long in Stan)
ts_data <- ts_data[seq(1, nrow(ts_data), dt/delta), ]

# Plotting
# Grid for x axis
x_grid <- seq(min(ts_data$state), max(ts_data$state), length.out = 500)

# ground truth cusp drift
ground_truth_drift <- cusp_drift(x_grid, r, alpha, beta, lambda, epsilon) %>% 
                      data.frame(x = x_grid, y = .)

# ground truth cusp diffusion
ground_truth_diff <- rep(0.5*cusp_dispersion(epsilon)^2, length(x_grid)) %>%
                     data.frame(x = x_grid, y = .)

# ground truth cusp density
ground_truth_stat_dens <- cusp_density(x_grid, r, alpha, beta, lambda, epsilon) %>% 
                          data.frame(x = x_grid, y = .)

#ggplot() +
#  geom_line(data = ground_truth_stat_dens, aes(x, y)) +
#  geom_line(aes(x_grid, nonparametric_stationary_density(ground_truth_drift$y, ground_truth_diff$y, x_grid)))

# find roots of drift
roots <- polyroot(c(-alpha, -beta, 0, 1))
roots <- roots + lambda
# remove complex roots
roots <- Re(roots[Im(roots) != 0])
roots %<>% sort()

# ground truth tipping point
real_root <- roots[2]

# calculate exit time
ground_truth_et_l <- exit_time_left_edge(ground_truth_drift$x, 
                                         which.min(abs(ground_truth_drift$x - real_root)), 
                                         ground_truth_drift$y, 
                                         ground_truth_diff$y)
ground_truth_et_r <- exit_time_right_edge(ground_truth_drift$x, 
                                          which.min(abs(ground_truth_drift$x - real_root)), 
                                          ground_truth_drift$y, 
                                          ground_truth_diff$y)
 
ground_truth_et <- data.frame("x" = c(ground_truth_et_l[[1]], ground_truth_et_r[[1]]), 
                              "et" = c(ground_truth_et_l[[2]], ground_truth_et_r[[2]]))

# calculate weighted mean exit time
p_st_left  <- ground_truth_stat_dens[1:which.min(abs(ground_truth_drift$x - real_root)), 2]
p_st_right <- ground_truth_stat_dens[which.min(abs(ground_truth_drift$x - real_root)):nrow(ground_truth_drift), 2]

mean_exit_time_left  <- mean_exit_time(p_st_left, ground_truth_et_l[[2]], ground_truth_et_l[[1]])
mean_exit_time_right <- mean_exit_time(p_st_right, ground_truth_et_r[[2]], ground_truth_et_r[[1]])
