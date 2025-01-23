
# Make part of figure 1 and set up for video comparing state space coverage with short time series and long time series

##################################################################################################################

# basic set up finding ground truth tipping point and range of the stationary density ############################

# set seed
seed <- 163
set.seed(seed)

# simulate short time series

# Run cusp catastrophe process
# cusp parameters
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# find tipping point
roots <- polyroot(c(-alpha, -beta, 0, 1))
roots <- roots + lambda
# remove complex roots
roots <- Re(roots[Im(roots) != 0])
roots %<>% sort()

# ground truth tipping point
tipping_point <- roots[2]

# find stationary density over a long range to see where it would be appropriate to cut it off

# find bounds of the cusp process
# Grid for x axis
x_grid <- seq(-2, 10, length.out = 500)

# calculate stationary density
ground_truth_stat_dens <- cusp_density(x_grid, r, alpha, beta, lambda, epsilon) %>% 
  data.frame(x = x_grid, y = .)

# find range of the data based on arbitrary probability cutoff

# maximum likelihood point (tallest peak of stationary density)
max_like <- max(ground_truth_stat_dens$y)
# see where it drops below 1000th of the largest likelihood value
#ground_truth_stat_dens$y < max_like / 1000
# find the indices of the range of the data where the stationary density is larger than max_like / 1000
bounds <- c(min(which((ground_truth_stat_dens$y < max_like / 1000) == F)), max(which((ground_truth_stat_dens$y < max_like / 1000) == F)))
# find length of the data range
d <- seq(ground_truth_stat_dens$x[bounds[1]], ground_truth_stat_dens$x[bounds[2]], length.out = 500)
data_range <- seq(d[1], d[length(d)], length.out = 500)

# calculate stationary density over the newly defined "full" data range
ground_truth_stat_dens <- cusp_density(data_range, r, alpha, beta, lambda, epsilon) %>% 
  data.frame(x = data_range, y = .)

##################################################################################################################

# calculate Kullback-Leibler divergence
# remember the ground truth is P and the comparison distribution is Q
KL_div <- function(P, Q) {
  
  return( sum(P*(log(P) - log(Q))) )
}

##################################################################################################################

# simulate 50 long time series of longest_sim_time total simulation time ########################################
num_time_series <- 50
longest_sim_time <- 500

df_lts = data.frame(matrix(vector(), longest_sim_time*2 + 1, num_time_series + 1))

# time step used for time series plot
dt <- 1

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- longest_sim_time*dt   # number of integer points
time       <- seq(from = 0, to = end_time, by = delta)
num_points <- 2*longest_sim_time+1

path <- here("output", "figures", "main_text", "figure 1")
dir.create(path)

if (file.exists(paste(path, "df_lts.rds", sep = "/"))) {

  df_lts <- readRDS(df_lts, file = paste(path, "df_lts.rds", sep = "/"))
} else {
  
  for (i in 1:num_time_series) {
    # Simulate from the cusp model
    lts_data <- cusp_euler_maruyama(y0 = NULL,
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
    df_lts[, i+1] <- lts_data[seq(1, nrow(lts_data), dt/(2*delta)), ]$state
  }
  
  df_lts[, 1] <- time[seq(1, nrow(lts_data), dt/(2*delta))]
  
  colnames(df_lts)[1] <- "time"
  
  # save
  saveRDS(df_lts, file = paste(path, "df_lts.rds", sep = "/"))
}

##################################################################################################################

# simulate 50 short time series of longest_sim_time total simulation time ########################################

df_sts = data.frame(matrix(vector(), longest_sim_time*2, num_time_series + 1))

# set up for the generating multiple short time series from the same process
time <- seq(from = 0, to = dt, by = delta)
# total number of points across all time series
total_points <- longest_sim_time*2
# number of short time series 
n_subjects <- longest_sim_time

df_sts[, 1] <- rep(time[seq(from = 1, to = length(time), by = floor(dt/delta))], longest_sim_time)

if (file.exists(paste(path, "df_sts.rds", sep = "/"))) {

  df_sts <- readRDS(df_sts, file = paste(path, "df_sts.rds", sep = "/"))
} else {
  
  for (i in 1:num_time_series) {
    
    ts_data <- data.frame(matrix(nrow = 0, ncol = 3)) 
    for (j in 1:n_subjects) {
      # Simulate from the cusp model
      short_ts_data <- cusp_euler_maruyama(y0 = NULL,
                                           times = time,
                                           seed = NULL,
                                           r = r,
                                           alpha = alpha,
                                           beta = beta,
                                           lambda = lambda,
                                           epsilon = epsilon) %>%
        cbind(time = time, subject = j) %>%
        as.data.frame() %>% 
        set_colnames(c("state", "time", "subject"))
      
      # subset data (otherwise it will take too long in Stan)
      ts_data <- rbind(ts_data, short_ts_data[seq(from = 1, to = length(time), by = floor(dt/delta)), ])
    }
    
    df_sts[, i+1] <- ts_data$state
  }
  
  colnames(df_sts)[1] <- "time"
  
  saveRDS(df_sts, file = paste(path, "df_sts.rds", sep = "/"))
}

##################################################################################################################

# now run analysis over all 50 time series ######################################################################

error_lts <- c()
error_sts <- c()

# calculate density errors over all 100 long time series
for (i in 1:longest_sim_time) {
  
  # find where to cut off time series for a given amount of time 
  time_bound <- min(which(as.integer(df_lts[,1]) == i))
  
  temp_error <- c()
  for (j in 1:num_time_series) {
    
    # density calculation
    approx_dens_lts <- density(df_lts[1:time_bound, j+1], 
                               n = nrow(ground_truth_stat_dens), 
                               from = data_range[1], 
                               to = data_range[length(data_range)])
    
    # visualise density comparison
    #ggplot() +
    #  geom_line(data = ground_truth_stat_dens, aes(x = x, y = y), size = 1) +
    #  geom_line(aes(x = approx_dens_lts$x, y = approx_dens_lts$y)) +
    #  theme_classic()
    
    # add a small number to any zeroes in the distribution in order to calculate the KL divergence
    approx_dens_lts$y[which(approx_dens_lts$y == 0)] <- approx_dens_lts$y[which(approx_dens_lts$y == 0)] + 
                                                        .Machine$double.eps
  
    # long time series error
    temp_error <- c(temp_error, KL_div(ground_truth_stat_dens$y, approx_dens_lts$y))
  }
  error_lts[[i]] <- temp_error
}
  
# calculate density errors over all 50 sets of short time series
# calculate density errors over all 50 long time series
for (i in 1:longest_sim_time) {
  
  # find where to cut off time series for a given amount of time 
  time_bound <- grep(1, df_sts[,1])[i]
  
  temp_error <- c()
  for (j in 1:num_time_series) {
    
    # density calculation
    approx_dens_sts <- density(df_sts[1:time_bound, j+1], 
                               n = nrow(ground_truth_stat_dens), 
                               from = data_range[1], 
                               to = data_range[length(data_range)])
    
    # visualise density comparison
    #ggplot() +
    #  geom_line(data = ground_truth_stat_dens, aes(x = x, y = y), size = 1) +
    #  geom_line(aes(x = approx_dens_sts$x, y = approx_dens_sts$y)) +
    #  theme_classic()
    
    # add a small number to any zeroes in the distribution in order to calculate the KL divergence
    approx_dens_sts$y[which(approx_dens_sts$y == 0)] <- approx_dens_sts$y[which(approx_dens_sts$y == 0)] + 
                                                        .Machine$double.eps
    
    # short time series error
    temp_error <- c(temp_error, KL_div(ground_truth_stat_dens$y, approx_dens_sts$y))
  }
  
  error_sts[[i]] <- temp_error
}

# get means and errors
lts_err_summary <- data.frame(matrix(vector(), longest_sim_time, 2))
sts_err_summary <- data.frame(matrix(vector(), longest_sim_time, 2))

for (i in 1:longest_sim_time) {
  
  lts_err_summary[i, 1] <- mean(error_lts[[i]]) 
  lts_err_summary[i, 2] <- sd(error_lts[[i]]) 
  
  sts_err_summary[i, 1] <- mean(error_sts[[i]]) 
  sts_err_summary[i, 2] <- sd(error_sts[[i]]) 
}

##################################################################################################################

# Plot and save ##################################################################################################

# normalise data
alpha_lts <- 1 / max(lts_err_summary[[1]])
alpha_sts <- 1 / max(sts_err_summary[[1]])

lts_err_summary[[1]] <- alpha_lts * lts_err_summary[[1]]
sts_err_summary[[1]] <- alpha_sts * sts_err_summary[[1]]
# make it increase towards 1 instead of decrease towards 0
lts_err_summary[[1]] <- 1 - lts_err_summary[[1]]
sts_err_summary[[1]] <- 1 - sts_err_summary[[1]]

# propagate error
lts_err_summary[[2]] <- alpha_lts * lts_err_summary[[2]]
sts_err_summary[[2]] <- alpha_sts * sts_err_summary[[2]]

# make dataframe for labelling plot
error_data <- data.frame("time" = rep(1:longest_sim_time, 2), 
                         "ts" = c(lts_err_summary[[1]], sts_err_summary[[1]]), 
                         "type" = c(rep("Long time series", longest_sim_time), rep("Short time series", longest_sim_time)),
                         "std_dev" = c(lts_err_summary[[2]], sts_err_summary[[2]]))

error_data$type <- factor(error_data$type, levels = c("Short time series", "Long time series"))

# for plotting with error
vals_lts <- which(error_data$type == "Long time series")
vals_sts <- which(error_data$type == "Short time series")

error_plot <- ggplot() +
  geom_line(data = error_data, aes(x = time, y = ts, colour = type, size = type, linetype = type)) +
  geom_ribbon(aes(x = error_data$time[vals_lts], 
                  ymin = error_data$ts[vals_lts] - error_data$std_dev[vals_lts]/sqrt(num_time_series), 
                  ymax = error_data$ts[vals_lts] + error_data$std_dev[vals_lts]/sqrt(num_time_series)), 
              fill = "gray", alpha = 0.4) + 
  geom_ribbon(aes(x = error_data$time[vals_sts], 
                  ymin = error_data$ts[vals_sts] - error_data$std_dev[vals_sts]/sqrt(num_time_series), 
                  ymax = error_data$ts[vals_sts] + error_data$std_dev[vals_sts]/sqrt(num_time_series)), 
              fill = "slateblue", alpha = 0.4) + 
  scale_colour_manual(values = c("slateblue", "gray")) +
  scale_size_manual(values = c(0.7, 0.9)) +
  scale_linetype_manual(values = c("solid", "solid")) +
  labs(x = "Simulation Time", y = "", title = "Agreement between data histogram and true stationary density") +
  theme_classic()

# create file to save output
path <- here("output", "figures", "main_text", "figure 1")
dir.create(path)
ggsave(error_plot, filename = paste(path, "short_vs_long_ts_plot.pdf", sep = "/"), dpi = 300, device = "pdf",
       height = 45, width = 90, units = "mm")

########################################################################################################################

# make histograms at t = 100
lts <- df_lts[1:(2*100+1), c(1, 4)]
sts <- df_sts[1:(2*100), c(1, 2)]

# density calculation
approx_dens_lts <- density(lts[, 2], 
                           n = nrow(ground_truth_stat_dens), 
                           from = data_range[1], 
                           to = data_range[length(data_range)])

approx_dens_sts <- density(sts[, 2], 
                           n = nrow(ground_truth_stat_dens), 
                           from = data_range[1], 
                           to = data_range[length(data_range)])

density_plot <- ggplot() +
  geom_line(data = ground_truth_stat_dens, aes(x, y), size = 2) + 
  geom_line(aes(x = approx_dens_sts$x, y = approx_dens_sts$y), colour = "slateblue", size = 1.2) +
  geom_line(aes(x = approx_dens_lts$x, y = approx_dens_lts$y), colour = "grey", size = 1.2) +
  theme_classic() +
  #scale_fill_manual(values = c("#D55E00", "#CC79A7")) +
  ylim(c(0, 0.8)) +
  xlim(c(0, 8)) +
  theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())

ggsave(density_plot, filename = paste(path, "densities.pdf", sep = "/"), dpi = 300, device = "pdf",
       height = 45, width = 45, units = "mm")



