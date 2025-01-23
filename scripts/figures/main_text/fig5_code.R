
# code to make fifth figure in the publication
########################################################################################################################
# run cusp model from figure 2
# code to make second figure in the publication

# Run cusp catastrophe process
# cusp parameters
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# set seed
seed <- 167 #231 #change seed back to 348757
set.seed(seed)

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- 1000   # number of integer points
time       <- seq(from = 0, to = end_time - 1, by = delta)
# time step used for time series plot
dt         <- 1
num_points <- ceiling(length(time) / (dt/delta))

# run cusp model
source(here("scripts", "cusp", "run_cusp.R"))

# set up for the generating multiple short time series from the same process
# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01
# time step used for time series plot 
dt         <- 0.2 #0.03*overall_et
# number of points per short time series
num_points <- 5
if (dt * num_points >= 1) {
  end_time <- ceiling(dt * num_points)
} else {
  end_time <- round(dt * num_points, digits = 2)
}
time <- seq(from = 0, to = end_time, by = delta)
# total number of points across all time series
total_points <- 250
# number of short time series 
n_subjects <- ceiling(total_points / num_points) # change 40 back to 80

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
  
  # subset data (otherwise it will take too long in Stan)
  ts_data <- rbind(ts_data, short_ts_data[seq(from = 1, to = length(time) - floor(dt/delta), by = floor(dt/delta)), ])
}

ground_truth <- 1
pred_inf <- 1
chosen_taxon <- 0

# create short time series plot
source(here("scripts", "model", "time_series_plot.R"))
ts_plot_ordered

# process data
pred_inf <- 1
ground_truth <- 1
ts_data_full <- ts_data
source(here("scripts", "model", "prepare_data.R"))

samples <- readRDS(here("output", "figures", "main_text", "figure 3", "fig3.1_data.rds"))$samples
data_read_from_file <- 1
source(here("scripts", "model", "run_drift_diff.R"))

# create file to save output
path <- here("output", "figures", "main_text", "figure 5")
dir.create(path)

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))
source(here("scripts", "model", "exit_time.R"))

# calculate exit time for drift and diffusion posterior means
cusp_et_l <- exit_time_left_edge(x_pred, pos_roots, drift[,1], diffusion[,1])
cusp_et_r <- exit_time_right_edge(x_pred, pos_roots, drift[,1], diffusion[,1])

cusp_et_tp_draws <- pos_exit_time_tp_CI
cusp_et_CI <- ET_centred_pos

# ground truth exit time calculation
# ground truth cusp drift
ground_truth_drift <- cusp_drift(x_pred, r, alpha, beta, lambda, epsilon) %>% 
  data.frame(x = x_pred, y = .)

# ground truth cusp diffusion
ground_truth_diff <- rep(0.5*cusp_dispersion(epsilon)^2, length(x_pred)) %>%
  data.frame(x = x_pred, y = .)

# ground truth cusp density
ground_truth_stat_dens <- cusp_density(x_pred, r, alpha, beta, lambda, epsilon) %>% 
  data.frame(x = x_pred, y = .)

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

cusp_et_plot <- ggplot() +
  geom_ribbon(aes(x = cusp_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = cusp_et_CI[, 1], ymax = cusp_et_CI[, 3]), fill = "mediumpurple2", alpha = 0.5) +
  geom_ribbon(aes(x = cusp_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = cusp_et_CI[, 1], ymax = cusp_et_CI[, 2]), fill = "mediumpurple3", alpha = 0.5) +
  geom_line(data = ground_truth_et, aes(x, et), size = 0.7, linetype = "dashed") +
  labs(title = "Cusp model", x = "System state", y = "Exit time") +
  geom_line(aes(c(cusp_et_l[[1]], cusp_et_r[[1]]), 
                c(cusp_et_l[[2]], cusp_et_r[[2]])), colour = "mediumpurple", size = 0.7) +
  #coord_cartesian(ylim = c(0, 180)) +
  theme_classic()

########################################################################################################################
# run exit time analysis for Prevotella

chosen_taxa <- "Prevotella_melaninogenica_et_rel"

ground_truth <- 0
pred_inf <- 1

# load data from figure 5
samples <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))[["samples"]]
x_pred  <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))[["x_pred"]]

data_read_from_file <- 1

# run analysis
source(here("scripts", "model", "run_drift_diff.R"))
source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))
source(here("scripts", "model", "exit_time.R"))

# calculate exit time for drift and diffusion posterior means
prevotella_et_l <- exit_time_left_edge(x_pred, pos_roots, drift[,1], diffusion[,1])
prevotella_et_r <- exit_time_right_edge(x_pred, pos_roots, drift[,1], diffusion[,1])

prevotella_et_tp_draws <- pos_exit_time_tp_CI
prevotella_et_CI <- ET_centred_pos

prevotella_et_plot <- ggplot() +
  geom_ribbon(aes(x = prevotella_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = prevotella_et_CI[, 1], ymax = prevotella_et_CI[, 3]), fill = "mediumpurple2", alpha = 0.5) +
  geom_ribbon(aes(x = prevotella_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = prevotella_et_CI[, 1], ymax = prevotella_et_CI[, 2]), fill = "mediumpurple3", alpha = 0.5) +
  labs(title = "Prevotella melaninogenica et rel.", x = "System state", y = "Exit time (days)") +
  geom_line(aes(c(prevotella_et_l[[1]], prevotella_et_r[[1]]), 
                c(prevotella_et_l[[2]], prevotella_et_r[[2]])), colour = "mediumpurple", size = 0.7) +
  coord_cartesian(ylim = c(0, 15000)) +
  theme_classic()


# plot tipping point region with drift and modality plot
figure <- ggarrange(cusp_et_plot + theme(text = element_text(size = 7)),
                    prevotella_et_plot + theme(text = element_text(size = 7)), 
                    ncol = 1, nrow = 2, align = "v")

# add common x-axis label
#figure <- annotate_figure(figure, bottom = textGrob("System state"))
# add title
#figure <- annotate_figure(figure, top = text_grob("Fig. 6: Exit time for real and simulated data", face = "bold", size = 14))
figure

# create file to save output
path <- here("output", "figures", "main_text", "figure 5")
dir.create(path)
# save
ggsave(paste(path, "fig5.pdf", sep = "/"), figure, device = "pdf", height = 130, width = 88, units = "mm", dpi = 500)













# Don't run the code below if you want to reproduce figure 5

########################################################################################################################
# run exit time analysis for lake Mendota
ground_truth <- 0
pred_inf <- 1
chosen_taxon <- 0

# load data from figure 3
samples <- readRDS(here("output", "figures", "main_text", "figure 3", "fig3_data.rds"))[["samples"]]
x_pred <- readRDS(here("output", "figures", "main_text", "figure 3",  "fig3_data.rds"))[["x_pred"]]

data_read_from_file <- 1

# run analysis
source(here("scripts", "model", "run_drift_diff.R")) #need to fix this
source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))
source(here("scripts", "model", "exit_time.R"))

# calculate exit time for drift and diffusion posterior means
mendota_et_l <- exit_time_left_edge(x_pred, pos_roots, drift[,1], diffusion[,1])
mendota_et_r <- exit_time_right_edge(x_pred, pos_roots, drift[,1], diffusion[,1])

mendota_et_tp_draws <- pos_exit_time_tp_CI
mendota_et_CI <- ET_centred_pos

# read in exit time paper's drift and diffusion
et_paper_drift <- readRDS(here("data", "exit_time_paper", "m_drift_app.rds"))
et_paper_diff <- readRDS(here("data", "exit_time_paper", "m_diff_app.rds"))

# calculate exit time based on the exit time paper's drift and diffusion
et_paper_root <- drift_roots(et_paper_drift$drift, et_paper_drift$x)$root[which(drift_roots(et_paper_drift$drift, et_paper_drift$x)$sign == 1)]
et_paper_et_l <- exit_time_left_edge(et_paper_drift$x, 
                                     et_paper_root, 
                                     et_paper_drift$drift, 
                                     et_paper_diff$diff)
et_paper_et_r <- exit_time_right_edge(et_paper_drift$x, 
                                      et_paper_root, 
                                      et_paper_drift$drift, 
                                      et_paper_diff$diff)

et_paper_et <- data.frame("x" = c(et_paper_et_l[[1]], et_paper_et_r[[1]]), 
                          "et" = c(et_paper_et_l[[2]], et_paper_et_r[[2]]))

mendota_et_plot <- ggplot() +
  geom_ribbon(aes(x = mendota_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = mendota_et_CI[, 1], ymax = mendota_et_CI[, 3]), fill = "mediumpurple2", alpha = 0.5) +
  geom_ribbon(aes(x = mendota_et_tp_draws$x[1:(length(x_pred) + 1)], ymin = mendota_et_CI[, 1], ymax = mendota_et_CI[, 2]), fill = "mediumpurple3", alpha = 0.5) +
  labs(title = "Lake Mendota", x = "System state", y = "Exit time (minutes)") +
  geom_line(aes(c(mendota_et_l[[1]], mendota_et_r[[1]]), 
                c(mendota_et_l[[2]], mendota_et_r[[2]])), colour = "mediumpurple", size = 0.7) +
  geom_line(data = et_paper_et, aes(x, et), colour = "royalblue4", size = 0.9, linetype = "dashed") +
  geom_vline(xintercept = most_prob_tp, size = 0.8, linetype = "dotted") +
  coord_cartesian(ylim = c(0, 18000)) +
  theme_classic()
