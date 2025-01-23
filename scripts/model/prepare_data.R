
# prepare drift and diffusion data for stan sampler
#if (!exists("ts_data")) {
#  ts_data <- otu_data[, c(5, 2, 4)]
#}

# change in x and y
ts_data$dt    <- NA
ts_data$dx    <- NA
ts_data$drift <- NA
ts_data$diff  <- NA

# check if continuous time series or short time series
if ("subject" %in% colnames(ts_data)) {
  
  # calculate changes
  for (i in 1:(nrow(ts_data)-1)) {
    if (ts_data$subject[i+1] == ts_data$subject[i]) {
    ts_data$dt[i]    <- ts_data$time[i+1] - ts_data$time[i]
    ts_data$dx[i]    <- ts_data$state[i+1] - ts_data$state[i]
    ts_data$drift[i] <- ts_data$dx[i] / ts_data$dt[i]
    ts_data$diff[i]  <- ts_data$dx[i]^2 / (2 * ts_data$dt[i])   # 1/2sigma^2
    } else next
  }
} else {
  
  # calculate changes
  for (i in 1:(nrow(ts_data)-1)) {
    ts_data$dt[i]    <- ts_data$time[i+1] - ts_data$time[i]
    ts_data$dx[i]    <- ts_data$state[i+1] - ts_data$state[i]
    ts_data$drift[i] <- ts_data$dx[i] / ts_data$dt[i]
    ts_data$diff[i]  <- ts_data$dx[i]^2 / (2 * ts_data$dt[i])   # 1/2sigma^2
  }
}

# remove NAs
ts_data  <- na.omit(ts_data)
# remove infinite values
ts_data <- ts_data[is.finite(ts_data$drift), ]
# re-order
ts_data <- ts_data[order(ts_data$state), ]
# if using log transform, can't have zero values. Should probably make it some small positive number so we don't
# lose data/information when it is already so limited
if (sum(which(ts_data$diff == 0)) > 0) {
  ts_data <- ts_data[-which(ts_data$diff == 0), ] 
}

# Stan data
n_chains <- 4
n_iter   <- 2000
c <- min(ts_data$state) + (max(ts_data$state) - min(ts_data$state)) / 2  # centre of the data

# grid for predictions (I don't recommend going below 100 points because then the resolution causes problems for other
# parts of the code)
x_pred <- seq(min(ts_data$state) - 0.3, max(ts_data$state) + 0.3, length.out = 100)

# data to be passed to Stan
stan_data <- list(N_real = nrow(ts_data),
                  N_pred = length(x_pred),
                  x_real = ts_data$state,
                  x_pred = x_pred,
                  dx = ts_data$dx,
                  dt = ts_data$dt,
                  c = c,
                  pred_inf = pred_inf)
