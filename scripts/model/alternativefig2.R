
# Lake Mendota
# load data
bga_data <- read.table(file = here("data", "Arani2021", "BGA_stdlevel_2011.csv"), sep = ",", header = T)

# plot original data
lake_men_plot <- ggplot(data = bga_data) +
  geom_line(aes(Tstep, X.stdlevel), colour = "royalblue4", size = 0.3) +
  labs(x = "Time (day number in 2011)", y = "Phycocyanin Level (log-transformed)", title = "Lake Mendota") +
  coord_cartesian(ylim = c(min(bga_data$X.stdlevel), max(bga_data$X.stdlevel))) +
  scale_x_continuous(breaks = seq(160, 250, 10)) +
  theme(text = element_text(size = 7)) +
  theme_classic()
#theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))
lake_men_plot


################################################################################
# Lake Mendota subsets 

# take stationary density from the exit time paper to draw starting values for lake Mendota
# I used an online tool to find points from the drift and diffusions function for lake Mendota
# in Arani et al. 2021 and saved them to the files data/ext_HITChip/Arani2021/mendota_drift.csv and mendota_diff.csv
# Then, I fit these points with a Gaussian process and that is what is being read in here
et_paper_drift <- readRDS(here("data", "Arani2021", "m_drift_app.rds"))
et_paper_diff <- readRDS(here("data", "Arani2021", "m_diff_app.rds"))

names(et_paper_drift) <- c("x", "y")
names(et_paper_diff) <- c("x", "y")

# calculate stationary density
stat_dens_et <- nonparametric_stationary_density(et_paper_drift$y, 
                                                 et_paper_diff$y, 
                                                 et_paper_drift$x) 

# basic parameters
ground_truth <- 0
pred_inf <- 1
chosen_taxon <- 0

# number of subjects
n_subjects <- 40
# number of points
n_points <- 5
# number of points in total
n_total <- n_points*n_subjects

# time step (we use 10 because in the exit time paper they had a time lag of 10)
dt <- 10

set.seed(236)
# y0 value drawn from regularly separated intervals
position <- floor(seq(from = 1, to = length(bga_data$X.stdlevel) - n_points*dt, length.out = n_subjects))

ts_data <- data.frame(matrix(nrow = 0, ncol = n_points))
for (i in 1:n_subjects) {
  
  # Subset lake mendota
  short_ts_data <- bga_data[seq(position[i], position[i] + (n_points - 1)*dt, by = dt), c(2, 3)] %>%
    cbind(subject = rep(i, n_points)) %>%
    set_colnames(c("time", "state", "subject"))
  
  # subset data (otherwise it will take too long in Stan)
  ts_data <- rbind(ts_data, short_ts_data)
}  

# add short time series to lake mendota time series plot
lake_men_plot <- lake_men_plot + geom_line(data = ts_data, aes(time, state, group = subject), size = 1, colour = "orchid1")

# convert to minutes
ts_data$time <- ts_data$time * 24 * 60

source(here("scripts", "model", "time_series_plot.R"))
ts_plot_ordered <- ts_data %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "orchid1") +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ts_plot_ordered

ts_data_full <- ts_data
source(here("scripts", "model", "prepare_data.R"))

samples <- readRDS(here("output", "figures", "main_text", "figure 2", "fig2_data.rds"))$samples
data_read_from_file <- 1

# if data_read_from_file = 0, MCMC result takes 2 hours to obtain on my laptop 
source(here("scripts", "model", "run_drift_diff.R"))

source(here("scripts", "model", "stability.R"))
if (length(neg_roots) > 1) {
  source(here("scripts", "model", "tipping_point.R"))
}

# potential function based off drift alone
deterministic_potential_function <- function(drift, x) {
  
  # change in x 
  dx <- (x[length(x)] - x[1]) / length(x)
  
  return(-cumsum(drift) * dx)
}

# Returns a dataframe of the stationary density with confidence intervals
# x        = vector, grid of points where the drift and diffusion were calculated
# samples  = stan output containing posterior information about drift and diffusion
potential_from_stan_samples_GP <- function(x_pred, samples, pred_inf) {
  
  if (pred_inf) {
    drift_draws <- rstan::extract(samples)$f_pred
  } else {
    drift_draws <- rstan::extract(samples)$f
  }
  
  df <- data.frame("x" = x_pred)
  
  # calculate stationary density for each drift and diffusion draw
  for (i in 1:nrow(drift_draws)) {
    df[, i+1] <- deterministic_potential_function(drift_draws[i, ], x_pred)
  }
  
  mean <- c(); upper_97.5 <- c(); lower_2.5 <- c(); upper_75 <- c(); lower_25 <- c()
  
  for (i in 1:nrow(df)) {
    mean[i] <- rowMeans(df[i, 2:ncol(df)], na.rm = T)
    #err <- c(mean[i] - 1.96*sd(df[i, 2:ncol(df)])/sqrt(nrow(drift_draws)), 
    #         mean[i] - 0.674*sd(df[i, 2:ncol(df)])/sqrt(nrow(drift_draws)),
    #         mean[i] + 0.674*sd(df[i, 2:ncol(df)])/sqrt(nrow(drift_draws)),
    #         mean[i] + 1.96*sd(df[i, 2:ncol(df)])/sqrt(nrow(drift_draws)))
    err <- quantile(df[i, 2:ncol(df)], probs = c(0.10, 0.25, 0.75, 0.90), na.rm = T)
    lower_2.5[i]  <- err[1]
    lower_25[i]   <- err[2]
    upper_75[i]   <- err[3]
    upper_97.5[i] <- err[4]
  }
  
  potential_posterior <- data.frame("x" = df$x, 
                                    "mean" = mean, 
                                    "lower_2.5" = lower_2.5,
                                    "lower_25" = lower_25,
                                    "upper_75" = upper_75,
                                    "upper_97.5" = upper_97.5)
  
  return(potential_posterior)
}

pos_potential <- potential_from_stan_samples_GP(x_pred, samples, pred_inf)

pot_et <- deterministic_potential_function(et_paper_drift$y, et_paper_drift$x)

figure1 <- ggarrange(lake_men_plot + geom_hline(yintercept = most_prob_tp, linetype = "dashed", colour = "orange")  + theme(text = element_text(size = 7)),
                     NULL,
                     ts_plot_ordered + coord_flip(xlim = c(min(bga_data$X.stdlevel), max(bga_data$X.stdlevel))) +
                       theme(axis.ticks.x = element_blank(),
                             axis.title.y = element_blank()) +
                       geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") +
                       scale_y_discrete(labels = rep("", n_subjects)) + theme(text = element_text(size = 7)), 
                     ncol = 3, nrow = 1, widths = c(2, 0.2, 1))

potential_plot <- ggplot(data = pos_potential) +
  geom_ribbon(aes(x = x_pred, ymin = lower_2.5, ymax = upper_97.5), fill = "cyan", alpha = 0.7) + 
  geom_ribbon(aes(x = x_pred, ymin = lower_25, ymax = upper_75), fill = "blue", alpha = 0.4) +
  labs(title = "", x = "System state", y = "Potential") +
  geom_line(aes(x, mean), alpha = 1, size = 0.7, colour = "royalblue3", show.legend = FALSE) +
  geom_line(aes(et_paper_drift$x, pot_et), size = 0.9, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-0.2, 0.2)) +
  guides(alpha = FALSE) +
  theme_classic()

potential_plot <- insert_xaxis_grob(potential_plot, CI_plot, grid::unit(.06, "null"), position = "top")

figure2 <- ggarrange(potential_plot,
                     stability_plot + theme(text = element_text(size = 7)),
                     ncol = 2, nrow = 1, align = "hv")

figure <- ggarrange(figure1, figure2, ncol = 1, nrow = 2, align = "hv", heights = c(1, 1.5))
figure

path <- here("output")
# save
ggsave(paste(path, "fig2redo.pdf", sep = "/"), figure, device = "pdf", height = 88, width = 88, units = "mm", dpi = 500)

