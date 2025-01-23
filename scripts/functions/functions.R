
####################### OUP FUNCTIONS START #############################

# Function to compute the exponential covariance kernel
# When you use the exponential kernel in a Gaussian process you get the OUP. 
exp_cov <- function(time, theta, diffusion) {

  # Initialise matrix
  covariance <- matrix(NA, nrow = length(time), 
                       ncol = length(time))

  for(i in 1:length(time)) {
    for(j in 1:length(time)) {
      covariance[i, j] = diffusion^2/(2*theta)*exp(-theta*abs(time[i]-time[j]));
    }
  }

  return(covariance)
}

# The stationary density of a OUP is a normal distribution determined by the mean and alpha: pi(x) = N(mean, alpha),
# where alpha is the variance. 
# (note that in R the standard deviation is usually used, instead of variance, e.g. rnorm(1, mean, sd))
OUP_stat_density <- function(x_grid, mu = 0, theta, diffusion) {
  
  y <- dnorm(x_grid, mu, diffusion/sqrt(2*theta))
  
  data.frame(x = x_grid, y = y)
}

# Drift is easy to compute as well
OUP_drift <- function(x_grid, mu = 0, theta) {
  
  y <- theta * (mu - x_grid)
  
  data.frame(x = x_grid, y = y)
}

####################### CUSP FUNCTIONS START #############################

# Constant dispersion
cusp_dispersion <- function(epsilon) {
  sqrt(epsilon)
}

# Drift
cusp_drift <- function(x, r, alpha, beta, lambda, epsilon) {
  r*(alpha + beta*(x - lambda) - (x - lambda)^3)
}

# Stationary density
cusp_density <- function(x, r, alpha, beta, lambda, epsilon) {
  
  unnormalized <- function(x, A = alpha, B = beta, L = lambda, E = epsilon, R = r) {
    exp( 2 * (A*(x - L) + .5*B*(x - L)^2 - .25*(x - L)^4 ) / (E/R))
  } 
  
  M <- integrate(unnormalized, -Inf, Inf)$value
  
  unnormalized(x)/M
}

# cusp simulator
cusp_euler_maruyama <- function(y0 = NULL, # intital value, if NULL --> randomly generates
                                r, alpha, beta, lambda, epsilon,
                                times, seed = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  if(is.null(y0)) {
    y <- sample_from_cusp(r, alpha, beta, lambda, epsilon)
  } else {
    y <- y0
  }
  
  for(i in 2:length(times)) {
    y_prev <- y[i-1]
    delta_t <- times[i] - times[i-1]
    
    a <- y_prev + cusp_drift(y_prev, r, alpha, beta, lambda, epsilon)*delta_t + 
         cusp_dispersion(epsilon)*rnorm(1, 0, sqrt(delta_t))
    y <- c(y, a)
    y
  }
  
  return(y)
}

# sample from the cusp stationary density
# cusp parameters
# x = range that covers the distribution
# N is the number of desired samples
# seed = seed
sample_from_cusp <- function(r, alpha, beta, lambda, epsilon, x = seq(0, 8, length.out = 1000), N = 1, seed = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # density without normalisation constant 
  prop_to <- function(x, R, A, B, L, E) {
    exp( 2 * (A*(x - L) + .5*B*(x - L)^2 - .25*(x - L)^4 ) / (E/R))
  }
  
  # find normalisation constant
  M <- integrate(prop_to, R = r, A = alpha, B = beta, L = lambda, E = epsilon, -Inf, Inf)$value  
  
  # density with normalisation constant
  PI <- function(x, M, R, A, B, L, E) {
    (1/M)*prop_to(x, R, A, B, L, E)
  }
  
  # calculate grid approximation of CDF
  CDF <- c()
  for (i in x) {
    CDF <- c(CDF, integrate(PI, M = M, R = r, A = alpha, B = beta, L = lambda, E = epsilon, -Inf, i)$value)
  }
  
  # generate uniform random sample and take the inverse CDF of it to generate samples
  unif_sample <- runif(N)
  S <- c()
  for (i in 1:N) {
    S[i] <- x[which.min(abs(CDF - unif_sample[i]))]
  }
  
  return(S)
}

# Calculate true positive and false negative rates for multiple simulations of many short time series each of the cusp process
#
# r, alpha, beta, lambda, epsilon: cusp parameters
# dt: time step
# num points: number of points per short time series
# total points: number of points added up across short time series
# num_sim: number of simulations to run. Each simulation generates new short time series
# num modes: is it a bistable cusp process or a unistable one
calculate_cusp_rates <- function(r = 0.2, alpha = -0.5, beta = 6, lambda = 4, epsilon = 1.2, dt = 0.2, num_points = 5, total_points = 250, num_sim = 100, num_modes = 2, seed = 100) {

  pred_inf <- 1
  ground_truth <- 0
  chosen_taxon <- 0
  data_read_from_file <- 0
  
  true_positives <- 0
  
  set.seed(seed)
  
  delta <- 0.01
  
  if (dt * num_points >= 1) {
    end_time <- ceiling(dt * num_points)
  } else {
    end_time <- round(dt * num_points, digits = 2)
  }
  time <- seq(from = 0, to = end_time, by = delta)
  
  # number of short time series 
  n_subjects <- ceiling(total_points / num_points)
  
  j <- 1
  while (j <= num_sim) {
    
    ts_data <- data.frame(matrix(nrow = 0, ncol = 3)) 
    for (k in 1:n_subjects) {
      # Simulate from the cusp model
      short_ts_data <- cusp_euler_maruyama(y0 = NULL,
                                           times = time,
                                           seed = NULL,
                                           r = r,
                                           alpha = alpha,
                                           beta = beta,
                                           lambda = lambda,
                                           epsilon = epsilon) %>%
        cbind(time = time, subject = k) %>%
        as.data.frame() %>% 
        set_colnames(c("state", "time", "subject"))
      
      # subset data (otherwise it will take too long in Stan)
      ts_data <- rbind(ts_data, short_ts_data[seq(from = 1, to = length(time) - floor(dt/delta), by = floor(dt/delta)), ])
    }
    
    source(here("scripts", "model", "prepare_data.R"), local = T)
    # run model ####################
    samples <- sampling(EM_model, stan_data, chains = 4, iter = 2000, cores = 4, seed = 234)
  
    # get posterior mean
    drift <- summary(samples, pars = "f_pred")[[1]][, c(1, 4:8)]
  
    # get posterior draws
    drift_draws <- rstan::extract(samples, "f_pred")[[1]]
    drift_draws <- drift_draws %>% 
      t %>%
      as.data.frame() %>% 
      mutate(x = x_pred)
    # end of model #################
    source(here("scripts", "model", "modality.R"), local = T)
    
    if (n_modes == num_modes) {
      true_positives <- true_positives + 1
    }
    
    print(paste0(paste0("Process is ", (j/num_sim)*100), "% done."))
    
    j <- j + 1
  }
  
  print("Remember: the true positive and false negative values are relative to the number of modes given.")
  return(list("true positives" = true_positives, "false negatives" = num_sim - true_positives))
}

####################### OTHER FUNCTIONS START #############################

transform_diff <- function(samples, pred_inf) {
  if (pred_inf) {
    diff_samples <- exp(rstan::extract(samples)$g_pred)
  } else {
    diff_samples <- exp(rstan::extract(samples)$g)
  }
  
  mean <- c(); upper_97.5 <- c(); lower_2.5 <- c(); upper_75 <- c(); lower_25 <- c()
  
  for (i in 1:ncol(diff_samples)) {
    mean[i]       <- mean(diff_samples[, i])
    err           <- quantile(diff_samples[, i], probs = c(0.025, 0.975), na.rm = T)
    lower_2.5[i]  <- err[1]
    upper_97.5[i] <- err[2]
    err           <- quantile(diff_samples[, i], probs = c(0.25, 0.75), na.rm = T)
    lower_25[i]   <- err[1]
    upper_75[i]   <- err[2]   
  }
  
  return(data.frame("mean" = mean, "lower_2.5" = lower_2.5, "lower_25" = lower_25, 
                    "upper_75" = upper_75, "upper_97.5" = upper_97.5))  
}

transform_drift <- function(samples) {
  drift_samples <- rstan::extract(samples)$f_pred / 10
  
  mean <- c(); upper_97.5 <- c(); lower_2.5 <- c(); upper_75 <- c(); lower_25 <- c()
  
  for (i in 1:ncol(drift_samples)) {
    mean[i]       <- mean(drift_samples[, i])
    err           <- quantile(drift_samples[, i], probs = c(0.025, 0.975), na.rm = T)
    lower_2.5[i]  <- err[1]
    upper_97.5[i] <- err[2]
    err           <- quantile(drift_samples[, i], probs = c(0.25, 0.75), na.rm = T)
    lower_25[i]   <- err[1]
    upper_75[i]   <- err[2]   
  }
  
  return(data.frame("mean" = mean, "lower_2.5" = lower_2.5, "lower_25" = lower_25, 
                    "upper_75" = upper_75, "upper_97.5" = upper_97.5))  
}

# Need to transform data so that it is of order unity for the model to work.
# Otherwise, variances and length scale are too large. Easier to transform data
# than priors. y is simply the y-values of the data
transform_to_unity <- function(y) {
  # if a decimal number
  if (max(y) < 1) {
    # truncate to one significant digit
    val <- signif(max(abs(y)), digits = 1)
    
    # find number of decimal places
    n <- nchar(strsplit(as.character(val), ".", fixed=TRUE)[[1]][2])
    
    return(y * 10 ^ n)
  } else {                             # if a non-decimal number
    # truncate to one significant digit
    val <- signif(max(y), digits = 1)
    
    return(y * 10 ^ (-(nchar(val)-1)))
  }
}

# potential function based off drift alone
deterministic_potential_function <- function(drift, x) {
  
  # change in x 
  dx <- (x[length(x)] - x[1]) / length(x)
  
  return(-cumsum(drift) * dx)
}

# Iacus page 75 (Stochastic Process with YUIMA)
# implemented math tricks to help with overflow and underflow
# returns dataframe of stationary density values and their x-values
# Calculate stationary density given drift and diffusion
# drift = (vector[real])   drift function
# diff  = (vector[real>0]) diffusion function
# x     = (vector[real])   x-values at which the drift and diffusion are evaluated
nonparametric_stationary_density <- function(drift, diff, x) {
  
  # change in x 
  dx <- (x[length(x)] - x[1]) / length(x)
  
  # Log scale measure
  integrand <- drift / (2 * diff)
  log_s     <- -2 * cumsum(integrand) * dx
  
  # Speed measure
  log_m     <- - log(2*diff) - log_s
  log_m_red <- log_m[-which.max(log_m)]
  
  log_dens  <- log_m - max(log_m) - log(dx) - log(1 + sum(exp(log_m_red - max(log_m))))
  
  # Normalize
  stationary_density <- exp(log_dens)
  
  return(stationary_density)
}

# Returns a dataframe of the stationary density with confidence intervals
# x        = vector, grid of points where the drift and diffusion were calculated
# samples  = stan output containing posterior information about drift and diffusion
nonparametric_stationary_density_from_stan_samples_GP <- function(x_pred, samples, pred_inf) {
  
  if (pred_inf) {
    drift_draws <- rstan::extract(samples)$f_pred
    diff_draws  <- 0.5*exp(rstan::extract(samples)$g_pred)
  } else {
    drift_draws <- rstan::extract(samples)$f
    diff_draws  <- 0.5*exp(rstan::extract(samples)$g)
  }

  df <- data.frame("x" = x_pred)
  
  # calculate stationary density for each drift and diffusion draw
  for (i in 1:nrow(drift_draws)) {
    df[, i+1] <- nonparametric_stationary_density(drift_draws[i, ], 
                                                  diff_draws[i, ], 
                                                  x_pred)
  }
  
  mean <- c(); upper_97.5 <- c(); lower_2.5 <- c(); upper_75 <- c(); lower_25 <- c()
  
  for (i in 1:nrow(df)) {
    mean[i] <- rowMeans(df[i, 2:ncol(df)], na.rm = T)
    err <- quantile(df[i, 2:ncol(df)], probs = c(0.025, 0.25, 0.75, 0.975), na.rm = T)
    lower_2.5[i]  <- err[1]
    lower_25[i]   <- err[2]
    upper_75[i]   <- err[3]
    upper_97.5[i] <- err[4]
  }
  
  density_posterior <- data.frame("x" = df$x, 
                                  "mean" = mean, 
                                  "lower_2.5" = lower_2.5,
                                  "lower_25" = lower_25,
                                  "upper_75" = upper_75,
                                  "upper_97.5" = upper_97.5)
  
  return(density_posterior)
}

stat_dens_modified_diff <- function(x_pred, drift_samples, diff_samples) {
  
  drift_draws <- rstan::extract(drift_samples)$f_pred
  diff_draws  <- exp(rstan::extract(diff_samples)$f_pred)
  
  df <- data.frame("x" = x_pred)
  
  # calculate stationary density for each drift and diffusion draw
  for (i in 1:nrow(drift_draws)) {
    df[, i+1] <- nonparametric_stationary_density(drift_draws[i, ], 
                                                  diff_draws[i, ], 
                                                  x_pred)
  }
  
  mean <- c(); upper_2.5 <- c(); lower_2.5 <- c()
  
  for (i in 1:nrow(df)) {
    mean[i] <- rowMeans(df[i, 2:ncol(df)], na.rm = T)
    err <- quantile(df[i, 2:ncol(df)], probs = c(0.025, 0.975), na.rm = T)
    upper_2.5[i] <- err[1]
    lower_2.5[i] <- err[2]
  }
  
  density_posterior <- data.frame("x" = df$x, 
                                  "mean" = mean, 
                                  "lower_2.5" = lower_2.5,
                                  "upper_2.5" = upper_2.5)
  
  return(density_posterior)
}

# Exit time ****************************************************************************************

# For a second order differential equation of the form ∂^2T/∂x^2 + p(x)∂T/∂x + q(x)T = r(x)
# It will be solved in matrix form Ay = b using the finite difference method O(dx^2)
# ys and yps are the boundary conditions
# If y values are known, ys = c(val1, val2). If mixed BCs, you can either have
# ys = c(val1, NA), yps = c(NA, val2) or ys = c(NA, val1), yps = c(val2, NA)
solveBVP <- function(p, q, r, x, ys, yps) {
  
  dx <- x[2] - x[1]                                  
  N  <- length(x)                
  A  <- matrix(0, nrow = N, ncol = N)  # ODE matrix
  b  <- r                              # RHS of matrix eqn
  
  # start of matrix
  A[1, 1] <- -2/dx^2 + q[1]
  A[1, 2] <- 1/dx^2 + p[1]/(2*dx)
  
  # body of matrix
  for (i in 2:(N-1)) {
    A[i, i-1] <- 1/dx^2 - p[i]/(2*dx)         # lower diagonal
    A[i, i]   <- -2/dx^2 + q[i]               # diagonal
    A[i, i+1] <- 1/dx^2 + p[i]/(2*dx)         # upper diagonal
  } 
  
  # end of matrix
  A[N, N-1] <- 1/dx^2 - p[N]/(2*dx)
  A[N, N] <- -2/dx^2 + q[N]
  
  # 2-point BCs
  if (sum(is.na(ys)) == 0) { 
    
    # remove rows and columns for BCs
    A <- A[-c(1, N), ]
    A <- A[, -c(1, N)]
    b <- b[-c(1, N)]
    
    # add BCs to b
    b[1] <- b[1] - (1/dx^2 - p[2]/(2*dx)) * ys[1]
    b[N-2] <- b[N-2] - (1/dx^2 + p[N-1]/(2*dx)) * ys[2]
    
  }
  
  # type mixed BCs
  if (T %in% is.na(ys)) {
    
    # remove rows and columns for BCs
    # case where the derivative is known at left boundary and the function at right
    if (is.na(ys[2])) {
      
      A <- A[-1, ]
      A <- A[, -1]
      b <- b[-1]
      
      # add BC to b
      b[1] <- b[1] - (1/dx^2 - p[2]/(2*dx)) * ys[1]
      # add other BC
      A[N-1, N-2] <- 2/dx^2                                       # a_(N,N-1) + a_(N,N+1)
      b[N-1] <- b[N-1] - 2*dx * (1/dx^2 + p[N]/(2*dx)) * yps[2]
      
      # case where the function is known at left boundary and the derivative at right
    } else if (is.na(ys[1])) {
      
      A <- A[-N, ]
      A <- A[, -N]
      b <- b[-N]
      
      # add BC to b
      b[N-1] <- b[N-1] - (1/dx^2 + p[N-1]/(2*dx)) * ys[2]
      # add other BC
      A[1, 2] <- 2/dx^2                                           # a_(1,0) + a_(1,2)
      b[1] <- b[1] + 2*dx * (1/dx^2 - p[1]/(2*dx)) * yps[1]
      
    } else {
      stop("Error: Improper boundary values.")
    }
  }
  
  # solve matrix equation
  y <- solve(A, b, tol = 1e-22)
  
  # return with BVs
  if (sum(is.na(ys)) == 0) {
    return(c(ys[1], y, ys[2])) 
  } else if (is.na(ys[2])) {
    return(c(ys[1], y)) 
  } else {
    return(c(y, ys[2]))
  }
}

# returns x_basins and exit time for a left edge 
# x_pred: vector of size N
# pos_roots: vector, location of roots in x_pred, integers
# drift: vector of size N
# diff: vector of size N
# start_point: position in x_pred vector you want the exit time calculation to start at
exit_time_left_edge <- function(x_pred, pos_roots, drift, diff, start_point = 1) {
  
  x <- c(); p <- c(); q <- c(); r <- c(); exit_time <- c()
  
  # find first basin
  p <- (drift[start_point:pos_roots[1]])/(diff[start_point:pos_roots[1]])
  q <- rep(0, pos_roots[1])
  r <- (-1)/(diff[start_point:pos_roots[1]])
  x <- x_pred[start_point:pos_roots[1]]
  
  ys <- c(NA, 0); yps <- c(1e-5, NA)
  
  exit_time <- solveBVP(p, q, r, x, ys, yps)
  
  return(list(x, exit_time))
}

# returns x_basins and exit time for intermediate basins
exit_time_intermediate <- function(x_pred, pos_roots, drift, diff) {
  
  x <- list(); p <- list(); q <- list(); r <- list(); exit_time <- list()
  
  ys <- c(0, 0)
  
  # find intermediate basins
  if (length(pos_roots) > 1) {
    for (j in 1:(length(pos_roots)-1)) {
      p[[j]] <- (drift[pos_roots[j]:pos_roots[j+1]])/(diff[pos_roots[j]:pos_roots[j+1]])
      q[[j]] <- rep(0, pos_roots[j+1] - pos_roots[j] + 1)
      r[[j]] <- (-1)/(diff[pos_roots[j]:pos_roots[j+1]])
      x[[j]] <- x_pred[pos_roots[j]:pos_roots[j+1]]
      
      exit_time[[j]] <- solveBVP(p[[j]], q[[j]], r[[j]], x[[j]], ys)
    }
  }
  return(list(x, exit_time))
}

# returns x_basins and exit time for a right edge 
# x_pred: vector of size N
# pos_roots: vector, location of roots in x_pred, integers
# drift: vector of size N
# diff: vector of size N
# end_point: position in x_pred vector you want the exit time calculation to stop at
exit_time_right_edge <- function(x_pred, pos_roots, drift, diff, end_point = length(x_pred)) {
  
  x <- c(); p <- c(); q <- c(); r <- c(); exit_time <- c()
  
  # find last basin
  p <- (drift[pos_roots[length(pos_roots)]:end_point]) / (diff[pos_roots[length(pos_roots)]:end_point])
  q <- rep(0, end_point - length(pos_roots) + 1)
  r <- (-1)/(diff[pos_roots[length(pos_roots)]:end_point])
  x <- x_pred[pos_roots[length(pos_roots)]:end_point]
  
  ys <- c(0, NA); yps <- c(NA, 1e-5);
  
  exit_time <- solveBVP(p, q, r, x, ys, yps)
  
  return(list(x, exit_time))
}

# compute Simpson's composite rule
composite_simpson <- function(x, f) {
  dx <- x[2] - x[1]
  N <- length(x)
  
  int_val <- (dx / 3) * (f[1] + 4 * sum(f[seq.int(1, N, 2)]) + 2 * sum(f[seq.int(2, N, 2)]) + f[N])
  
  return(int_val)
}

# to find mean exit time, where the boundary conditions are the edges of the basin of attraction.
# Next, calculate T_av = ∫p_st(x) / ∫p_st(y)dy T(x) dx where both integrals run over the basin of attraction so a_∫^b 

# p_st is the stationary density
# t is the exit time
# x is range of the basin
mean_exit_time <- function(p_st, t, x) {
  
  c <- composite_simpson(x, p_st)    # weight for the stationary density
  f <- t * p_st / c                  # function to be integrated
  
  t_av <- composite_simpson(x, f)    # weighted-mean exit time
  
  return(t_av)
}

# make routine that finds roots of drift, differentiating them by positive or negative slope
drift_roots <- function(drift, x_pred) {
  
  d_roots <- data.frame("root" = NA, "sign" = NA)
  
  # find where it crosses the x-axis
  crossings <- round(zero_crossings(drift), digits = 0)
  
  # check that there are roots so that it doesn't return an error
  if (any(is.na(crossings))) {
    return(d_roots)
  } else {
    # the difference between indices can't be one
    if (1 %in% base::diff(crossings)) {
      #crossings[which(base::diff(crossings) == 1) + 1] <- crossings[which(base::diff(crossings) == 1) + 1] + 1
      crossings[which(base::diff(crossings) == 1)] <- crossings[which(base::diff(crossings) == 1)] - 1
    }
    
    for (i in unique(crossings)) {
      # make sure not an endpoint since can't calculate slope on both sides
      if (i != 1 && i != length(x_pred)) {
        # for determining whether it has a positive or negative slope
        slope_right <- (drift[i+1] - drift[i]) / (x_pred[i+1] - x_pred[i])
        slope_left <- (drift[i] - drift[i-1]) / (x_pred[i] - x_pred[i-1])
        
        # normal root
        if (sign(slope_right) == sign(slope_left)) {
          d_roots <- rbind(d_roots, data.frame("root" = i, "sign" = sign(slope_right)))
        } else if (sign(slope_right) != sign(slope_left) && drift[i] != 0) {
          # it thinks it's a bounce but really it should be the two points beside the "bounce"
          d_roots <- rbind(d_roots, data.frame("root" = c(i-1, i+1), "sign" = c(sign(slope_left), sign(slope_right))))
        } else {
          # an actual bounce
          d_roots <- rbind(d_roots, data.frame("root" = i, "sign" = 0))
        }
      } else if (i == 1) {
        # deal with left endpoint
        slope_right <- (drift[i+1] - drift[i]) / (x_pred[i+1] - x_pred[i])
        d_roots <- rbind(d_roots, data.frame("root" = i, "sign" = sign(slope_right)))
      } else if (i == length(x_pred)) {
        # deal with right endpoint
        slope_left <- (drift[i] - drift[i-1]) / (x_pred[i] - x_pred[i-1])
        d_roots <- rbind(d_roots, data.frame("root" = i, "sign" = sign(slope_left)))
      }
    }
    
    # remove first row of NAs
    d_roots <- d_roots[-1, ]
    
    return(d_roots)
  }
}

# HitChip functions

## Data getters ****************************** ####

# These functions read and prepare the data for gaussian process interpolators

# There are dublicated time indices in hitchip data. 
# This functions selects the sample with higher read count and returns the edited phyloseq
get_edit_hitchip_pseq <- function(pseq) {
  
  # Edit taxa_names
  atlas_taxa <- taxa_names(pseq) %>% 
    gsub(" ", "_", .)
  taxa_names(pseq) <- atlas_taxa
  
  # Add read counts
  sample_data(pseq) <- cbind(meta(pseq), read_count = rowSums(t(abundances(pseq))))
  
  # Subjects with multiple samples
  multiple_id <- meta(pseq) %>% 
    filter(time != 0) %>% 
    pull(subject) %>% 
    unique()
  
  # Mark all longitudinal samples as duplicates, some of these are removed
  sample_data(pseq) <- cbind(meta(pseq), dup = meta(pseq) %>% 
                               mutate(dup = ifelse(subject %in% multiple_id, "yes", "no")) %>% 
                               pull(dup))
  
  meta <- meta(pseq)
  
  for(i in unique(meta$subject)) {
    
    # Time indices
    subject_time <- meta %>% 
      filter(subject == i) %>% 
      pull(time)
    
    # Duplicated times  
    dup_time <- subject_time[which(duplicated(subject_time))]
    
    # Remove smaller read counts from meta
    for(t in unique(dup_time)) {
      max_read_count <- meta %>% 
        filter(subject == i, time == t) %>% 
        pull(read_count) %>% 
        max()
      
      meta[meta$subject == i &
             meta$time == t &
             (meta$read_count != max_read_count), ] <- NA
      
      meta <- meta %>% drop_na
    }
  }
  
  sample_data(pseq) <- meta
  
  return(pseq)
}

pick_baseline <- function (A) {
  
  # THIS function picks the baseline subset of complete Atlas for close inspection
  # print("Pick sample subset (adult fecal samples with RBB extraction with no reported health issues)")
  
  tmp <- subset(meta(A),
                age_group == "adult" &
                  sample_type == "fecal" &
                  dna_extraction_method == "rbb" &
                  !health_status == "compromised")
  
  # Pick a single baseline sample per subject
  # s <- rownames(tmp[!duplicated(tmp[, "subject"]),])
  # pseq.baseline <- prune_samples(s, A)
  A.baseline <- baseline(prune_samples(rownames(tmp), A))
  
  A.baseline
  
}

pick_timeseries <- function (A, atlas.metadata.full) {
  
  # THIS function picks the subset of complete Atlas that includes only subjects with time series from
  # accepted projects and criteria
  
  # Ignore some projects from time series analysis due to 
  # Include only the selected adult samples
  # poor quality or intervention issues
  
  # Pick time series
  # Use only Fecal RBB samples from non-compromised subjects
  # with no reported antibiotics or medication
  # Also remove specified projects due to questionable sample quality
  
  m0 <- subset_samples(A,
                       (dna_extraction_method == "rbb" &
                          sample_type == "fecal" &
                          # Also include compromised so we can possibly compare
                          #(!health_status == "compromised" | is.na(health_status)) &
                          (antibio == 0 | is.na(antibio)) &
                          (medication == 0 | is.na(medication)))
  )
  
  m <- sample_data(m0)
  
  # Some comments
  # "ProvasPlus", # Intervention but accepted bc we showed earlier no effect https://www.ncbi.nlm.nih.gov/pubmed/21829582
  #"IHO_IBS",             # 21 (anti-TNF medication but some untreated)
  #"SNAS Patients Study", # All compromised
  #"SYMPPIS",             # 137 x (probiotic); all compromised
  
  rmproj <- c(
    # Remove these as intervention studies
    "ELIPA", # n=100 Intervention / (weight loss)
    "IBSVEGAN",  # n = 38 (gut cleansing + ancient yoga diet)
    "PURGING", # Intervention
    "HEALTHGRAIN",         # n = 104; Intervention (grains); all compromised
    "TURKU PET STUDY",     # Intervention; Obese; mostly compromised
    
    # Filtered above already but listed here if criteria change
    # as these are not good quality studies for time series analysis
    "Extraction",
    "POUCH", 
    "PROTEOMICS",
    "Spatial",
    "STORAGE PROJECT"
  )
  m <- subset(m, !project %in% rmproj)
  
  # Identify subjects with multiple time points
  s <- names(which(lapply(split(m$time, m$subject), function (x) {length(unique(na.omit(x)))}) > 1))
  
  # Accepted samples for longitudinal analyses
  # Includes both adults and babies
  s <- rownames(subset(m, subject %in% s))
  
  # Summary table of the available time series
  ts <- atlas.metadata.full[s, ] %>% group_by(project) %>% summarise(subjects = length(unique(subject)), samples = n())
  print(ts)
  
  # Prepare organized time series data
  A.timeseries <- prune_samples(s, A) # was: atlas.time.pseq
  
  return(A.timeseries)
  
}


#' @title Exporting phyloseq Data in CSV Files
#' @description Writes the otu, taxonomy and metadata in csv files.
#' @param x \code{\link{phyloseq-class}} object
#' @param type 'OTU' or 'TAXONOMY' or 'METADATA'
#' @param path Path to the directory/folder where the data will be written.
#' Uses the working directory by default.
#' @return  Output file path (a string)
#' @seealso read_phyloseq
#' @export
#' @examples 
#' #data(dietswap)
#' #pseq <- dietswap
#' ## By default writes all info at once (ie OTU/TAXONOMY/METADATA)
#' #write_phyloseq(pseq) 
#' #write_phyloseq(pseq, 'OTU')
#' #write_phyloseq(pseq, 'TAXONOMY')
#' #write_phyloseq(pseq, 'METADATA')
#' @keywords utilities
write_phyloseq <- function(x, type="all", path=getwd()) {
  
  .Deprecated("", "The microbiome::write_phyloseq will be deprecated in a future release.")
  
  type <- toupper(type)
  
  # TODO make read_phyloseq as well
  if (type == "OTU" || type == "all") {
    f <- paste(path, "otu_table.csv", sep="/")
    message("Writing OTU in the file ", f)
    # y <- as.data.frame(x@otu_table);
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    # Let us use abundances function here as it is guaranteed to be taxa x
    # samples always
    y <- abundances(x)
    write.csv(y, file=f, fileEncoding="UTF-16LE")
    
  }
  
  if (type == "TAXONOMY" || type == "all") {
    # Renamed from TAXA to TAXONOMY as the latter is used elsewhere
    f <- paste(path, "taxonomy_table.csv", sep="/")
    message("Writing TAXONOMY in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- as.data.frame(tax_table(x))
    write.csv(y, file=f, fileEncoding="UTF-16LE")
    
  }
  
  if (type == "METADATA" || type == "all") {
    f <- paste(path, "metadata_table.csv", sep="/")
    message("Writing METADATA in the file ", f)
    if (f %in% dir(path)) {
      warning(paste("The file with the same name", f,
                    "exists in the given path and is overwritten."))
    }
    y <- meta(x)
    write.csv(y, file=f, fileEncoding="UTF-16LE")
  }
  
  return(path)
  
}



select_first_timepoint <- function (A) {
  
  # Subjects with time information
  M1 <- filter(meta(A), !is.na(subject) & !is.na(time)) %>%
    arrange(time) 
  # Find first occurrence (time point) for each subject
  M1 <- M1[match(unique(M1$subject), M1$subject), ]
  
  # Subjects with no time information: exclude those with multiple samples
  M2 <- filter(meta(A), !is.na(subject) & is.na(time)) %>%
    group_by(subject) %>%
    mutate(n = n()) %>%
    filter(n == 1)
  
  # Accepted unique first time point samples
  samples <- c(M1$sample, M2$sample)
  
  subset_samples(A, sample %in% samples)
  
}

filter_atlas <- function (keep, atlas.metadata = NULL, atlas.sampleinfo = NULL, atlas) {
  
  if (!is.null(atlas.metadata)) {
    atlas.metadata <- atlas.metadata[keep,]
  }
  
  if (!is.null(atlas.sampleinfo)) {  
    atlas.sampleinfo <- atlas.sampleinfo[keep,]
  }
  
  atlas$oligo <- atlas$oligo[, keep]
  for (lev in c("species", "L0", "L1", "L2")) {
    for (meth in names(atlas[[lev]])) {
      atlas[[lev]][[meth]] <- atlas[[lev]][[meth]][, keep]
      if (lev == "L2") {
        rownames(atlas[[lev]][[meth]]) <- gsub("^Clostridium\ \\(sensu stricto\\)$", "Clostridium_sensu_stricto", rownames(atlas[[lev]][[meth]]))
      }
      # Remove periods
      rownames(atlas[[lev]][[meth]]) <- gsub("_$", "", gsub("_+", "_", gsub(" ", "_", gsub("\\.", "_", rownames(atlas[[lev]][[meth]])))))
      # Put back 'rel.' periods
      rownames(atlas[[lev]][[meth]]) <- gsub("rel_$", "rel", rownames(atlas[[lev]][[meth]]))
    }
  }
  for (lev in names(atlas)[grep("remap", names(atlas))]) {
    atlas[[lev]] <- atlas[[lev]][, keep]
  }
  
  list(atlas.metadata = atlas.metadata, atlas.sampleinfo = atlas.sampleinfo, atlas = atlas)
  
}




# Add age group info for the samples
age.cohort <- function (em.info) {
  
  age.group <- rep(NA, nrow(em.info))
  names(age.group) <- rownames(em.info)
  
  age.group[which(em.info$age >= 0 & em.info$age < 1)] <- "age_0"
  age.group[which(em.info$age >= 1 & em.info$age < 2)] <- "age_1"
  age.group[which(em.info$age >= 2 & em.info$age < 3)] <- "age_2"
  age.group[which(em.info$age >= 3 & em.info$age < 6)] <- "age_3_to_5"
  age.group[which(em.info$age >= 6 & em.info$age < 7)] <- "age_6"
  age.group[which(em.info$age >= 7 & em.info$age < 8)] <- "age_7"
  age.group[which(em.info$age >= 8 & em.info$age < 13)] <- "age_8_to_12"
  age.group[which(em.info$age >= 13 & em.info$age < 16)] <- "age_13_to_15"
  age.group[which(em.info$age >= 16 & em.info$age < 18)] <- "age_16_to_17"
  age.group[which(em.info$age >= 18 & em.info$age < 20)] <- "age_18_to_19"  
  age.group[which(em.info$age >= 20 & em.info$age < 30)] <- "age_20s"
  age.group[which(em.info$age >= 30 & em.info$age < 40)] <- "age_30s"
  age.group[which(em.info$age >= 40 & em.info$age < 50)] <- "age_40s"
  age.group[which(em.info$age >= 50 & em.info$age < 60)] <- "age_50s"
  age.group[which(em.info$age >= 60 & em.info$age < 70)] <- "age_60s"
  age.group[which(em.info$age >= 70 & em.info$age < 80)] <- "age_70s"
  age.group[which(em.info$age >= 80)] <- "age_80_plus"
  
  age.group <- factor(age.group, levels = unique(age.group[order(em.info$age)]))
  
  age.group
  
}




random_pairs <- function (x) {
  
  # Generate random sample pairs
  
  # Pick metadata
  m <- meta(x)
  
  # Random pairing indices
  m$pairingID <- sample(rep(1:(nrow(m)/2), 2))[1:nrow(m)]
  
  # Pick random pairs
  keep <- names(which(table(m$pairingID) > 1))
  m <- subset(m, pairingID %in% keep)
  split.random <- split(m$sample, m$pairingID)
  split.random <- lapply(split.random, function (x) {sample(x, 2)})
  
  # Remove samples from the same subject
  rminds <- which(sapply(split.random, function (x) {na.omit(m[x[[1]],"subject"] == m[x[[2]], "subject"])}) == TRUE)
  if (length(rminds) > 0) {
    split.random <- split.random[setdiff(names(split.random), names(rminds))]
  }
  
  split.random
  
}

within_subject_pairs <- function (x) {
  
  # Pick metadata
  m <- meta(x)
  m$pairingID <- m$subject
  
  # Pick two random samples per subject
  # TODO consider also temporal differences
  keep <- names(which(table(m$pairingID) > 1))
  m <- subset(m, pairingID %in% keep)
  split.person <- split(m$sample, m$pairingID)
  split.person <- lapply(split.person, function (x) {sample(x, 2)})
  
  split.person
  
}




# Add age group info for the samples
age.cohort.numeric <- function (em.info) {
  
  n <- nrow(em.info)
  age.group <- rep(NA, n)
  names(age.group) <- rownames(em.info)
  
  age.group[which(em.info$age >= 0 & em.info$age < 1)]   <- 0
  age.group[which(em.info$age >= 1 & em.info$age < 2)]   <- 1
  age.group[which(em.info$age >= 2 & em.info$age < 3)]   <- 2
  age.group[which(em.info$age >= 3 & em.info$age < 4)]   <- 3
  age.group[which(em.info$age >= 4 & em.info$age < 5)]   <- 4
  age.group[which(em.info$age >= 5 & em.info$age < 6)]   <- 5
  age.group[which(em.info$age >= 6 & em.info$age < 7)]   <- 6
  age.group[which(em.info$age >= 7 & em.info$age < 8)]   <- 7
  age.group[which(em.info$age >= 8 & em.info$age < 10)]  <- 8
  age.group[which(em.info$age >= 10 & em.info$age < 16)] <- 10
  age.group[which(em.info$age >= 15 & em.info$age < 20)] <- 15
  age.group[which(em.info$age >= 20 & em.info$age < 30)] <- 20
  age.group[which(em.info$age >= 30 & em.info$age < 40)] <- 30
  age.group[which(em.info$age >= 40 & em.info$age < 50)] <- 40
  age.group[which(em.info$age >= 50 & em.info$age < 60)] <- 50
  age.group[which(em.info$age >= 60 & em.info$age < 70)] <- 60
  age.group[which(em.info$age >= 70 & em.info$age < 80)] <- 70
  age.group[which(em.info$age >= 80)] <- 80
  
  age.group
  
}



# Add age group info for the samples
bmi.class <- function (em.info) {
  
  bmi.group <- rep(NA, nrow(em.info))
  names(bmi.group) <- rownames(em.info)
  
  bmi <- em.info$bmi
  
  # Split obese in three groups
  bmi.group <- rep(NA, nrow(em.info))
  bmi.group[bmi < 18.5] <- "underweight"
  bmi.group[bmi >= 18.5 & bmi < 25] <- "lean"
  bmi.group[bmi >= 25] <- "overweight"
  bmi.group[bmi >= 30] <- "obese"
  bmi.group[bmi >= 35] <- "severe"
  bmi.group[bmi >= 40] <- "morbid"
  bmi.group[bmi >= 45] <- "super"
  
  bmi.group <- factor(bmi.group, levels = c("underweight", "lean", "overweight", "obese", "severe", "morbid", "super"))
  
  # Ensure BMI group only assigned to adults
  keep <- setdiff(1:length(bmi.group), which(em.info$age >= 18))
  bmi.group[keep] <- NA
  
  bmi.group
  
}



plot_atlas_profile <- function (pseq, x, y) {
  
  df <- data.frame(sample_data(pseq)[, x])
  df$xvar <- df[[x]]  
  
  if (y %in% names(sample_data(pseq))) {
    df$signal <- sample_data(pseq)[[y]]
  } else if (y %in% rownames(otu_table(pseq))) {
    df$signal <- as.vector(otu_table(pseq)[y,])
  }
  
  # Randomize sample order to avoid visualization biases
  df <- df[sample(nrow(df)),]
  
  # Order the factor levels (x variables) by their standard deviation
  df$xvar <- factor(df$xvar, levels = names(sort(sapply(split(df$signal, df$xvar), sd))))
  
  # Order by x variable
  # (randomization affects the ordering within each factor level if this is a factor)
  df <- df[order(df$xvar),]
  
  df$index <- 1:nrow(df)
  
  library(ggplot2)
  theme_set(theme_bw(20))
  p <- ggplot(df, aes(x = index, y = signal, color = xvar))
  p <- p + geom_point()
  p <- p + guides(color = guide_legend(ncol = 2, title = x))
  p <- p + scale_y_log10()
  p <- p + ggtitle(y)
  p <- p + xlab("Sample index")
  p <- p + ylab(y)    
  
  p
  
}



pick_cross_sectional <- function (meta, subject.field = "unique.subjectID", time.field = "timepoint") {
  # Use the first time point if time information is given;
  # otherwise pick one of the samples at random
  set.seed(3424)
  s <- split(meta, meta[[subject.field]])
  ids <- sapply(s, function (x) {
    ind <- order(as.numeric(x[[time.field]]))[[1]];
    if (all(is.na(x$time))) {ind <- sample(nrow(x), 1)}; 
    rownames(x)[[ind]]
  })
  meta <- meta[ids,]
  rownames(meta)
}

select_samples <- function (meta, remove.compromised = TRUE) {
  
  # print("Pick fecal RBB samples")
  library(dplyr)
  meta <- filter(meta, sample_type == "fecal" & dna_extraction_method == "rbb")
  rownames(meta) <- meta$anonymized.sampleID
  
  # Remove compromised samples
  if (remove.compromised) {
    compromised <- rownames(subset(meta, health_status %in% "compromised"))
    meta <- meta[setdiff(rownames(meta), compromised), ]
  }
  
  # print("Remove outlier samples with strikingly increased oligo median")
  # as in our Nat Comm Paper
  outliers <- names(which(apply(atlas.full$oligo, 2, median) > 1.9))
  keep <- setdiff(rownames(meta), outliers)
  meta <- meta[keep, ]
  
  # Remove subjects with antibiotics
  meta <- meta[-which(meta$antibio == 1), ]
  
  # Remove FMT samples (all healthy samples in TURN TRIAL are FMT; see groupcodes)
  meta <- filter(meta, !project == "TURN TRIAL")
  
  # Remove high BMI
  # meta <- meta[-which(meta$bmi > 43), ]
  
  # print("Order subjects")
  meta <- meta[order(meta$unique.subjectID),]
  
  # print("Pick the corresponding HITChip Atlas data and Metadata")
  selected.samples <- meta$anonymized.sampleID
  
  # data.date <- date()
  # print("Write data")
  # write.csv(meta, file = "Atlas1000.tab", quote = FALSE, row.names = FALSE)
  # save(atlas, meta, phylogeny.info.filtered, phylogeny.info.full, data.date, file = paste("Atlas1000-", gsub(" ", "-", date()), ".RData", sep = ""), compress = "xz")
  
  selected.samples
  
}


subject_time_shifts <- function (df, subject = "unique.subjectID", sample = "anonymized.sampleID", time = "time", signal = "signal") {
  
  # Pick data subset to analyze
  dfsub <- df[, c(subject, sample, time, signal)]
  colnames(dfsub) <- c("subject", "sample", "time", "signal")
  
  # Ensure ordering by time
  library(dplyr)
  dfsub <- arrange(dfsub, time)
  
  # Pick data for each subject separately
  spl <- split(dfsub, dfsub$subject)
  
  tabs <- NULL
  
  for (subj in names(spl)) {
    
    # Pick time and signal
    times <- as.numeric(spl[[subj]]$time)
    signal <- as.numeric(spl[[subj]]$signal)
    
    # Remove possible duplicates
    keep <- !duplicated(times)
    times <- times[keep]
    signal <- signal[keep]
    
    # List of combinations of t1 vs. t2
    # where t2 > t1
    n <- length(times)
    grid <- expand.grid(1:n, 1:n)
    grid <- grid[grid[, 2] > grid[, 1],]
    if (length(grid) == 0) {return(NULL)}
    colnames(grid) <- c("ind1", "ind2")
    
    shifts <- NULL
    for (k in 1:nrow(grid)) { 
      i <- grid[k, 1]
      j <- grid[k, 2]
      tp1 <- times[[i]]
      tp2 <- times[[j]]
      s1 <- signal[[i]]
      s2 <- signal[[j]]
      # step: denotes difference between time points
      # (1 can be later used to pick just the consecutive time point shifts)
      x <- c(t1 = tp1, t2 = tp2, s1 = s1, s2 = s2, step = j-i)
      shifts <- rbind(shifts, x)
    }
    
    library(dplyr)
    shifts <- as.data.frame(shifts)
    shifts <- mutate(shifts, ds = s2 - s1)
    shifts <- mutate(shifts, dt = t2 - t1)
    shifts$subject <- subj
    
    tabs <- rbind(tabs, shifts)
  }
  
  tabs
  
}



#' Polish time series data frame
#'
#' @param meta data.frame with the following fields: subject, sample, time, and optionally other columns
#' @param min.timepoints Include only subjects with at least this many timepoints
#' @param subject.field Name of the subjectID field
#' @param time.field Name of the time field
#' @param sample.field Name of the sampleID field
#' @return Polished data.frame where the first time points in each 
#' 	  subject are shifted such that the first time point is zero. 
#'	  The samples are ordered in temporal order within each subject. 
#'        Samples with no time information are removed
polish_timeseries <- function (meta, min.timepoints = 1, subject.field = "unique.subjectID", time.field = "time", sample.field = "anonymized.sampleID") {
  
  # Remove samples with no time data		  
  meta <- meta[!is.na(meta[[time.field]]), ]
  
  # List subjects who have time series
  spl <- split(meta[[time.field]], meta[[subject.field]])
  time.subjects <- names(which(sapply(spl, function (x) {length(unique(x))>1})))
  
  # Pick samples from subjects who have real time series
  df <- meta[, c(time.field, subject.field, sample.field)]
  names(df) <- c("time", "subject", "sample")
  df <- subset(df, subject %in% time.subjects)
  
  spl <- split(df, df$subject)
  df2 <- NULL
  for (subj in names(spl)) {
    spl2 <- spl[[subj]]
    times <- as.numeric(as.character(spl2$time))
    
    # Shift the times such that the first timepoint is zero
    spl2$time <- times - min(times)
    
    # Order samples within subject based on time
    o <- order(spl2$time)    	      
    spl2 <- spl2[o, ]
    
    # If there are multiple baseline points pick the first one and ignore the rest
    inds <- which(spl2$time == 0)
    if (length(inds) > 1) {
      spl2 <- spl2[-inds[-1],]
    }
    
    df2 <- rbind(df2, spl2)
  }
  
  # Keep subjects with sufficiently many non-zero timepoints
  keep.subjects <- names(which(sort(table(df2[df2$time>0,]$subject)) >= min.timepoints))
  
  df2[df2$subject %in% keep.subjects, ]
  df2 <- df2[, c("sample", "subject", "time")]
  rownames(df2) <- df2$sample
  # Return the original column names?
  #names(df2) <- c(sample.field, subject.field, time.field)
  
  df2
  
}

group.nationalities <- function (nationality.orig) {
  
  nationality2 <- as.character(nationality.orig) # rep(NA, length(nationality.orig))
  nationality2[nationality.orig %in% c("AAM", "USA")] <- "USA"
  nationality2[nationality.orig %in% c("FIN", "SWE", "NOR")] <- "Scandinavia"
  nationality2[nationality.orig %in% c("GBR", "IRL")] <- "UKIE"
  nationality2[nationality.orig %in% c("ESP", "ITA", "FRA")] <- "SouthEur"
  nationality2[nationality.orig %in% c("DK", "DEU", "BEL", "NLD")] <- "CentralEur"
  nationality2[nationality.orig %in% c("POL", "SRB", "CC")] <- "EastEur"
  nationality2[nationality.orig %in% c("PAK", "JAP", "IND", "IDN")] <- "Asia" 
  
  factor(nationality2)
  
}




# Estimate the number of modes for each row of a data matrix
# mydat: phylotypes vs. samples
nmodes <- function (x, det.th = 5, bw.adjust = 1, min.density = 30) {
  x <- na.omit(unname(x))
  a <- livpotential_ews(x, bw = "nrd", detection = det.th, min.density = min.density)
  a$max.points
}


