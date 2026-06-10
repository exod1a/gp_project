
# code to make supplementary figure for revisions

# first we have to calculate t_c
# In order to do that, first we need to calculate the distance in the formula t_c = d / <v>
# calculate characteristic time scale based on model parameters
r       <- 0.1
alpha   <- 0
beta    <- 6
lambda  <- 4
epsilon <- 0.9

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

# result is t_c = 75.59745
t_c <- (data_range^2) / epsilon

# now we have to generate short time series
# set seed
seed <- 13
set.seed(seed)

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- 100   # number of integer points
time       <- seq(from = 0, to = end_time - 1, by = delta)
# time step used for time series plot
dt         <- 1
num_points <- ceiling(length(time) / (dt/delta))

# run cusp model
source(here("scripts", "cusp", "run_cusp.R"))

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
# to scale delta, the resolution of the simulation
target_k <- 40
aim_delta <- 0.01
# time step used for time series plot 
dt         <- 0.01*t_c
# set max cutofff
delta <- min(dt / target_k, aim_delta)
k <- round(dt / delta)
delta <- dt / k
# number of points per short time series
num_points <- 5
if (dt * num_points >= 1) {
  end_time <- ceiling(dt * num_points)
} else {
  end_time <- round(dt * num_points, digits = 2)
}
time <- seq(from = 0, to = end_time, by = delta)
# total number of points across all time series
total_points <- 200
# number of short time series 
n_subjects <- ceiling(total_points / num_points) 

ts_data <- data.frame(matrix(nrow = 0, ncol = 3)) 
for (i in 1:n_subjects) {
  # Simulate from the cusp model
  short_ts_data <- cusp_euler_maruyama(y0 = NULL,
                                       times = time,
                                       seed = NULL,
                                       r = r,
                                       alpha = alpha,
                                       beta = beta,
                                       lambda = lambda,
                                       epsilon = epsilon) %>%
    cbind(time = time, subject = i) %>%
    as.data.frame() %>% 
    set_colnames(c("state", "time", "subject"))
  
  short_ts_data <- short_ts_data[seq(from = 1, to = length(time), by = floor(dt/delta)), ]
  short_ts_data <- short_ts_data[1:num_points, ]
  # subset data (otherwise it will take too long in Stan)
  ts_data <- rbind(ts_data, short_ts_data)
}

# create short time series plot
source(here("scripts", "model", "time_series_plot.R"))

# look at transitions
ts_data <- ts_data %>%
  arrange(subject, time) %>%
  group_by(subject) %>%
  mutate(
    prev_state = lag(state),
    
    # Detect crossings
    cross_up   = prev_state < real_root & state >= real_root,
    cross_down = prev_state >= real_root & state < real_root
  ) %>%
  summarise(
    transition = case_when(
      any(cross_up, na.rm = TRUE) & !any(cross_down, na.rm = TRUE) ~ "left_to_right",
      any(cross_down, na.rm = TRUE) & !any(cross_up, na.rm = TRUE) ~ "right_to_left",
      any(cross_up, na.rm = TRUE) & any(cross_down, na.rm = TRUE)  ~ "both",
      TRUE ~ "none"
    ),
    .groups = "drop"
  ) %>%
  right_join(ts_data, by = "subject")

ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val), colour = transition), size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ts_hist <- ggplot() + 
  geom_histogram(data = ts_data, aes(x=state), colour="black", fill="#f0d4ec", bins = 22) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "State", y = "Count") +
  theme_classic()
ts_hist

# create file to save output
path <- here("output", "figures", "extended_data", "figure 4")
dir.create(path)

# run model
# process data
pred_inf <- 1
ground_truth <- 1
ts_data_full <- ts_data
source(here("scripts", "model", "prepare_data.R"))

if (file.exists(here("output", "figures", "extended_data", "figure 4", "fig4.1_data.rds"))) {
  
  samples <- readRDS(here("output", "figures", "extended_data", "figure 4", "fig4.1_data.rds"))$samples
  data_read_from_file <- 1
  source(here("scripts", "model", "run_drift_diff.R"))
} else {
  
  data_read_from_file <- 0
  source(here("scripts", "model", "run_drift_diff.R"))
  
  fig4_data <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(fig4_data, file = paste(path, "fig4.1_data.rds", sep = "/"))
}

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))

drift_plot <- drift_plot + geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-2, 2)) + theme(text = element_text(size = 7))

drift_plot <- insert_xaxis_grob(drift_plot, CI_plot, grid::unit(.06, "null"), position = "top")

ts_data_full$transition <- factor(ts_data_full$transition,
                        levels = c("left_to_right", "right_to_left", "none"))

final_ts_plot <- ts_data_full %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val), colour = transition), size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  scale_color_manual(values = c("orchid1", "orchid1", "royalblue4")) +
  theme_classic() + 
  geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

# run analyses on Puhti over 100 simulations with these cusp parameters, 200 points, 0.01t_c, unbiased
unbiased_stats <- read.table(here("output", "figures", "extended_data", "figure 4", "fig4.1_data_puhti.txt"), header = T)
n_bis_ub <- unbiased_stats$n_bis
unbiased_stats <- unbiased_stats[, -5]

# run analyses on Puhti over 100 simulations with these cusp parameters, ~200 points (some removed for biasing), 0.01t_c, biased
# for some reason, one simulation failed
biased_stats <- read.table(here("output", "figures", "extended_data", "figure 4", "fig4.2_data_puhti.txt"), header = T)
n_bis_b <- biased_stats$n_bis
biased_stats <- biased_stats[, -5]

# convert to long format for box plot
unbiased_stats_long <- unbiased_stats %>%
                       pivot_longer(cols = everything(),
                                    names_to = "variable",
                                    values_to = "value") %>%
                       mutate(group = "Unbiased")

# convert to long format for box plot
biased_stats_long <- biased_stats %>%
                     pivot_longer(cols = everything(),
                                  names_to = "variable",
                                  values_to = "value") %>%
                     mutate(group = "Biased")

# Combine
stats_long <- bind_rows(unbiased_stats_long, biased_stats_long)
stats_long$group <- factor(stats_long$group, levels = c("Unbiased", "Biased"))

# look for significant differences
# add star in Inkscape
pvals <- stats_long %>%
  group_by(variable) %>%
  summarise(p = t.test(value ~ group)$p.value) %>%
  mutate(star = case_when(
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  ))

# dist_right is one-star significant. Add this to plot

# Boxplot with significance
stats_boxplots <- ggplot(stats_long, aes(x = variable, y = value, fill = group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  theme_classic() +
  labs(x = "", y = "Error", fill = "Coverage") #+
  #scale_fill_manual(values = c("Unbiased" = "#1f78b4", "Biased" = "#33a02c"))

# plot histogram of proportion of draws that were bistable over the simulations
n_bis_ub <- as.data.frame(n_bis_ub) %>% mutate(group = "Unbiased")
colnames(n_bis_ub) <- c("n_bis", "group")
n_bis_b <- as.data.frame(n_bis_b) %>% mutate(group = "Biased")
colnames(n_bis_b) <- c("n_bis", "group")
n_bis <- bind_rows(n_bis_ub, n_bis_b)
n_bis$group <- factor(n_bis$group, levels = c("Unbiased", "Biased"))

# Overlapping histograms with transparency
n_bis_plot <- ggplot(n_bis, aes(x = n_bis, fill = group)) +
  geom_histogram(position = "identity",   # overlay instead of stacking
                 alpha = 0.5,            # transparency (0 = fully transparent, 1 = opaque)
                 bins = 20,              # number of bins
                 color = "black") +      # optional: add border around bins
  #scale_fill_manual(values = c("Unbiased" = "#1f78b4", "Biased" = "#33a02c")) +
  theme_classic() +
  labs(x = "Percentage of Bistable Posterior Draws", y = "Count", fill = "Coverage")

# plot tipping point region with drift and modality plot
figure1 <- plot_grid(final_ts_plot, 
                     drift_plot, 
                     n_bis_plot,
                     ncol = 1, nrow = 3, align = "hv")

# the code to make the biased plots: transitions only from left to right
# now we have to generate short time series
# set seed
seed <- 17
set.seed(seed)

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- 100   # number of integer points
time       <- seq(from = 0, to = end_time - 1, by = delta)
# time step used for time series plot
dt         <- 1
num_points <- ceiling(length(time) / (dt/delta))

# run cusp model
source(here("scripts", "cusp", "run_cusp.R"))

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
# to scale delta, the resolution of the simulation
target_k <- 40
aim_delta <- 0.01
# time step used for time series plot 
dt         <- 0.01*t_c
# set max cutofff
delta <- min(dt / target_k, aim_delta)
k <- round(dt / delta)
delta <- dt / k
# number of points per short time series
num_points <- 5
if (dt * num_points >= 1) {
  end_time <- ceiling(dt * num_points)
} else {
  end_time <- round(dt * num_points, digits = 2)
}
time <- seq(from = 0, to = end_time, by = delta)
# total number of points across all time series
total_points <- 200
# number of short time series 
n_subjects <- ceiling(total_points / num_points) 

ts_data <- data.frame(matrix(nrow = 0, ncol = 3)) 
for (i in 1:n_subjects) {
  # Simulate from the cusp model
  short_ts_data <- cusp_euler_maruyama(y0 = NULL,
                                       times = time,
                                       seed = NULL,
                                       r = r,
                                       alpha = alpha,
                                       beta = beta,
                                       lambda = lambda,
                                       epsilon = epsilon) %>%
    cbind(time = time, subject = i) %>%
    as.data.frame() %>% 
    set_colnames(c("state", "time", "subject"))
  
  short_ts_data <- short_ts_data[seq(from = 1, to = length(time), by = floor(dt/delta)), ]
  short_ts_data <- short_ts_data[1:num_points, ]
  # subset data (otherwise it will take too long in Stan)
  ts_data <- rbind(ts_data, short_ts_data)
}

# create short time series plot
source(here("scripts", "model", "time_series_plot.R"))

# look at transitions
ts_data <- ts_data %>%
  arrange(subject, time) %>%
  group_by(subject) %>%
  mutate(
    prev_state = lag(state),
    
    # Detect crossings
    cross_up   = prev_state < real_root & state >= real_root,
    cross_down = prev_state >= real_root & state < real_root
  ) %>%
  summarise(
    transition = case_when(
      any(cross_up, na.rm = TRUE) & !any(cross_down, na.rm = TRUE) ~ "left_to_right",
      any(cross_down, na.rm = TRUE) & !any(cross_up, na.rm = TRUE) ~ "right_to_left",
      any(cross_up, na.rm = TRUE) & any(cross_down, na.rm = TRUE)  ~ "both",
      TRUE ~ "none"
    ),
    .groups = "drop"
  ) %>%
  right_join(ts_data, by = "subject")

ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val), colour = transition), size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

ts_hist <- ggplot() + 
  geom_histogram(data = ts_data, aes(x=state), colour="black", fill="#f0d4ec", bins = 22) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "State", y = "Count") +
  theme_classic()
ts_hist

# remove short time series that transition from right to left and go both ways
# This way, it is biased: only transitions from left to right
to_remove <- which(ts_data$transition %in% c("both", "right_to_left"))
ts_data <- ts_data[-to_remove, ]

ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val), colour = transition), size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# run model
# process data
pred_inf <- 1
ground_truth <- 1
ts_data_full <- ts_data
source(here("scripts", "model", "prepare_data.R"))

if (file.exists(here("output", "figures", "extended_data", "figure 4", "fig4.2_data.rds"))) {
  
  samples <- readRDS(here("output", "figures", "extended_data", "figure 4", "fig4.2_data.rds"))$samples
  data_read_from_file <- 1
  source(here("scripts", "model", "run_drift_diff.R"))
} else {
  
  data_read_from_file <- 0
  source(here("scripts", "model", "run_drift_diff.R"))
  
  fig6_data <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(fig6_data, file = paste(path, "fig4.2_data.rds", sep = "/"))
}

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))

drift_plot <- drift_plot + geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-2, 2)) + theme(text = element_text(size = 7))

drift_plot <- insert_xaxis_grob(drift_plot, CI_plot, grid::unit(.06, "null"), position = "top")

final_ts_plot <- ts_data_full %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val), colour = transition), size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  scale_color_manual(values = c("orchid1", "royalblue4")) +
  theme_classic() + 
  geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

# plot tipping point region with drift and modality plot
figure2 <- plot_grid(final_ts_plot, 
                     drift_plot, 
                     stats_boxplots, 
                     ncol = 1, nrow = 3, align = "hv")

# full figure
figure <- plot_grid(figure1, 
                    figure2,
                    ncol = 2, nrow = 1, align = "hv")

# save
ggsave(paste(path, "fig4.pdf", sep = "/"), figure, device = "pdf", height = 180, width = 200, units = "mm", dpi = 500)

