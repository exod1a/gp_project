
# calculate the drift function like they did in the exit time paper
et_drift_calc <- function(state_values, bins, N_lags, dt) {
  D1 <- c()
  state_app <- c()
  
  for (j in 1:length(bins)) {
    res <- drift_lags(state_values, bins, j, N_lags)
    M1_data <- data.frame("dt" = (1:N_lags)*dt, "M1" = res[[1]])
    # fit regression for bin j 
    if (N_lags > 1) {
      D1 <- c(D1, coef(lm(data = M1_data, M1 ~ dt))[2])
    } else {
      D1 <- c(D1, (M1_data$M1) / (M1_data$dt))
    }
    state_app <- c(state_app, res[[2]])
  }
  
  return(list(D1, state_app))
}

# calculate N_lag number of drift values at a given been to be fit with a linear regression
drift_lags <- function(state_values, bins, x_value, N_lags) {
  M1_app <- c()
  
  state_app <- mean(state_values[bins[[x_value]]])
  
  # The points for the different time lags
  for (i in 1:N_lags) {
    M1_app <- c(M1_app, sum(state_values[bins[[x_value]] + i] - state_values[bins[[x_value]]]) / (length(bins[[x_value]])))
  }
  
  return(list(M1_app, state_app))
}

# Lake Mendota data

# binning
ts_data <- bga_data
ts_data <- ts_data[, -1]
names(ts_data) <- c("time", "state")

bin_size <- 1500
dt <- (ts_data$time[2] - ts_data$time[1]) * 24 * 60  # convert from days to minutes 
N_lags <- 10

# find points that fall in each bin
bin_x_vals <- list()
# remove last ten points from binning so that there will always be N_lag time points left over to 
# calculate the drifts
to_remove <- c((nrow(ts_data)-(N_lags - 1))):(nrow(ts_data))

for (i in 1:(floor(nrow(ts_data) / bin_size))) {
  bin_x_vals[[i]] <- which(ts_data$state[-to_remove] %in% sort(ts_data$state[-to_remove])[((i-1)*bin_size+1):(i*bin_size)])
}

state_values <- ts_data$state; bins <- bin_x_vals
res <- et_drift_calc(ts_data$state, bin_x_vals, N_lags, dt)

D1 <- res[[1]]
state_app <- res[[2]]

# compare to their drift result
mendota_drift <- read.table(file = here("data", "Arani2021", "mendota_drift.csv"), sep = ",", header = T)

ggplot() +
  geom_point(aes(state_app, D1), size = 2, alpha = 0.7, shape = 1) +
  geom_line(data = mendota_drift, aes(x, y), colour = "red", size = 0.7) +
  geom_hline(yintercept = 0) +
  theme_classic()
  

########################################################################################################################

# get ts_data from fig. 2

bin_size <- 1
dt <- (ts_data$time[2] - ts_data$time[1]) * 24 * 60  # convert from days to minutes 
N_lags <- 1

# find points that fall in each bin
bin_x_vals <- list()
# remove last ten points from binning so that there will always be N_lag time points left over to 
# calculate the drifts
to_remove <- c((nrow(ts_data)-(N_lags - 1))):(nrow(ts_data))

for (i in 1:(floor(nrow(ts_data) / bin_size))) {
  bin_x_vals[[i]] <- which(ts_data$state[-to_remove] %in% sort(ts_data$state[-to_remove])[((i-1)*bin_size+1):(i*bin_size)])
}

state_values <- ts_data$state; bins <- bin_x_vals
res <- et_drift_calc(ts_data$state, bin_x_vals, N_lags, dt)

D1 <- res[[1]]
state_app <- res[[2]]

# compare to their drift result
mendota_drift <- read.table(file = here("data", "Arani2021", "mendota_drift.csv"), sep = ",", header = T)

ggplot() +
  geom_point(aes(state_app, D1), size = 2, alpha = 0.7, shape = 1) +
  geom_line(data = mendota_drift, aes(x, y), colour = "red", size = 0.7) +
  geom_hline(yintercept = 0) +
  theme_classic()











# Cusp model data

# run cusp model first to get ts_data
# binning

bin_size <- 100
dt <- (ts_data$time[2] - ts_data$time[1]) 
N_lags <- 10

# find points that fall in each bin
bin_x_vals <- list()
# remove last ten points from binning so that there will always be N_lag time points left over to 
# calculate the drifts
to_remove <- c((nrow(ts_data)-(N_lags - 1))):(nrow(ts_data))

for (i in 1:(floor(nrow(ts_data) / bin_size))) {
  bin_x_vals[[i]] <- which(ts_data$state[-to_remove] %in% sort(ts_data$state[-to_remove])[((i-1)*bin_size+1):(i*bin_size)])
}

res <- et_drift_calc(ts_data$state, bin_x_vals, N_lags, dt)
res2 <- et_drift_calc(ts_data$state, bin_x_vals, 1, dt)
res3 <- et_drift_calc(ts_data$state, bin_x_vals, 2, dt)

D1 <- res[[1]]
state_app <- res[[2]]

D1_ours <- res2[[1]]
state_app_ours <- res2[[2]]

D1_other <- res3[[1]]
state_app_other <- res3[[2]]

ggplot() +
  geom_point(aes(state_app, D1), size = 2, alpha = 0.7, shape = 1) +
  geom_point(aes(state_app_ours, D1_ours), size = 2, alpha = 0.7, shape = 1, colour = "red") +
  geom_point(aes(state_app_other, D1_other), size = 2, alpha = 0.7, shape = 1, colour = "purple") +
  geom_line(data = c_drift, aes(x, y), size = 0.7) + labs(title = "", x = "State", y = "Drift") + 
  geom_hline(yintercept=0, size = 0.4) +
  geom_point(aes(roots, 0), size = 2.5) + 
  coord_cartesian(xlim = c(-2, 2), ylim = c(-1, 1)) +
  theme_classic()
