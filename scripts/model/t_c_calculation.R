
# calculate the characteristic time scale for the cusp model and the bimodal unistable model

############################################################################################################
# calculate characteristic time scale using many long cusp realisations
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# set seed
seed <- 231 
set.seed(seed)

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- 200   # number of integer points
time       <- seq(from = 0, to = end_time - 1, by = delta)
# time step used for time series plot
dt         <- 1
num_points <- ceiling(length(time) / (dt/delta))

char_ts <- c()

for (j in 1:100) {
  
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
  
  
  ts_data_full <- ts_data
  source(here("scripts", "model", "prepare_data.R"))
  
  char_ts[j] <- (max(ts_data_full$state) - min(ts_data_full$state)) / mean(2*ts_data$diff)
  print(j)
}

hist(char_ts)
mean(char_ts)

############################################################################################################

# first we need to calculate the distance in the formula t_c = d / <v>
# calculate characteristic time scale based on model parameters
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# Grid for x axis
x_grid <- seq(-2, 10, length.out = 500)

# calculate stationary density
ground_truth_stat_dens <- cusp_density(x_grid, r, alpha, beta, lambda, epsilon) %>% 
  data.frame(x = x_grid, y = .)

# plot
ggplot(data = ground_truth_stat_dens) +
  labs(title = "", x = "State", y = "Density") +
  geom_ribbon(aes(x, ymin = 0, ymax = y), fill = "#f0d4ec", alpha = 1, colour = "black") +
  theme_classic()

# maximum likelihood point (tallest peak of stationary density)
max_like <- max(ground_truth_stat_dens$y)
# see where it drops below 1000th of the largest likelihood value
ground_truth_stat_dens$y < max_like / 1000
# find the indices of the range of the data where the stationary density is larger than max_like / 1000
bounds <- c(min(which((ground_truth_stat_dens$y < max_like / 1000) == F)), max(which((ground_truth_stat_dens$y < max_like / 1000) == F)))
# find length of the data range
d <- seq(ground_truth_stat_dens$x[bounds[1]], ground_truth_stat_dens$x[bounds[2]], length.out = 500)
data_range <- d[length(d)] - d[1]


##########################
# Ville's method

t_c <- (data_range^2) / epsilon

##########################################################################################################################################
# calculate t_c for bimodal case

sample_from_dens <- function(x, dens, N = 1, seed = NULL) {
  
  # set seed
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  # change in x 
  dx <- diff(x)[1]
  
  # calculate grid approximation of CDF
  CDF <- cumsum(dens)*dx
  
  # generate uniform random sample and take the inverse CDF of it to generate samples
  unif_sample <- runif(N)
  S <- c()
  for (i in 1:N) {
    S[i] <- x[which.min(abs(CDF - unif_sample[i]))]
  }
  
  # return sample(s)
  return(S)
}

# run custom stochastic process
run_SDE <- function(DRIFT, DIFFUSION, t, ground_truth_stat_dens, seed = NULL, y0 = NULL) {
  
  if(!is.null(seed)) {
    set.seed(seed)
  }
  
  dt <- diff(t)
  
  if(is.null(y0)) {
    y <- sample_from_dens(ground_truth_stat_dens[[1]], ground_truth_stat_dens[[2]])
  } else {
    y <- y0
  }
  
  for(i in 2:length(t)) {
    y_prev <- y[i-1]
    y_next <- y_prev + sapply(y_prev, DRIFT) * dt[i-1] + sqrt(sapply(y_prev, DIFFUSION)) * rnorm(1, 0, sqrt(dt[i-1]))
    y <- c(y, y_next)
  }
  
  return(y)
}

# custom drift function
DRIFT <- function(x) {
  
  FUN1 <- function(x) {
    exp(-0.08*x) - 0.95
  }
  FUN2 <- function(x) {
    -0.5*x^2 + 0.05
  }
  
  if (length(x) > 1) {
    x1 <- x[x <= 0]
    x2 <- x[x > 0]
    
    y1 <- sapply(x1, FUN1)
    y2 <- sapply(x2, FUN2)
    
    return(c(y1, y2))
  } else {
    if (x < 0) return(sapply(x, FUN1))
    if (x > 0) return(sapply(x, FUN2))
  }
}

# custom diffusion function
DIFFUSION <- function(x) {
  
  FUN1 <- function(x) {
    2*0.422*exp(-(2*x - 0.6)^2) #0.8*exp(-(x-1)^2) #0.5*exp(3*x)
  }
  
  FUN2 <- function(x) {
    2*0.422*exp(-(2*x - 0.6)^2) #0.1*sin(5*x) + 0.2943 #0.35*atan(4*x) + 0.5 #0.3*sin(5*x) + 0.5
  }
  
  
  if (length(x) > 1) {
    x1 <- x[x <= 0]
    x2 <- x[x > 0]
    
    y1 <- sapply(x1, FUN1)
    y2 <- sapply(x2, FUN2)
    
    return(c(y1, y2))
  } else {
    if (x < 0) return(sapply(x, FUN1))
    if (x > 0) return(sapply(x, FUN2))
  }
}


x_grid <- seq(-2, 2, length.out = 5000)

# custom density plot
ground_truth_stat_dens <- data.frame("x" = x_grid, 
                                     "y" = nonparametric_stationary_density(DRIFT(x_grid), 
                                                                            0.5*DIFFUSION(x_grid), 
                                                                            x_grid))

# maximum likelihood point (tallest peak of stationary density)
max_like <- max(ground_truth_stat_dens$y)
# see where it drops below 1000th of the largest likelihood value
ground_truth_stat_dens$y < max_like / 1000
# find the indices of the range of the data where the stationary density is larger than max_like / 1000
bounds <- c(min(which((ground_truth_stat_dens$y < max_like / 1000) == F)), max(which((ground_truth_stat_dens$y < max_like / 1000) == F)))
# find length of the data range
d <- seq(ground_truth_stat_dens$x[bounds[1]], ground_truth_stat_dens$x[bounds[2]], length.out = 500)
data_range <- d[length(d)] - d[1]

x_grid <- seq(d[1], d[length(d)], length.out = 5000)

t_c_bimodal <- (data_range^2) / (mean(DIFFUSION(x_grid)))
t_c_bimodal

# t_c used in calculations was 13.53338

########################################################################################################################


# calculate characteristic time scale for Prevotella
chosen_taxa <- "Prevotella_melaninogenica_et_rel"

# full time series
ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
names(ts_data_full)[4] <- "state"

data_range <- (max(ts_data_full$state) - min(ts_data_full$state))

# read in time series data with changes e.g. dx, dt
ts_data <- readRDS(here("data", "ext_HITChip", chosen_taxa, "taxon_data.rds"))

diffusion <- (ts_data$delta)^2 / (ts_data$dt)

t_c <- (data_range^2) / mean(diffusion)

# check if it is close to 45
0.01*t_c

# in reality 0.02t_c
