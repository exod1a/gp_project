
# code to make extended data figure 1 in the publication
# takes roughly 5 minutes to run

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
# read in data from points taken from their graph and fitted with a GP
# note that these are not really "ground truth" functions. They are just named that way
# to go with how I have set up the code
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
n_subjects <- 10
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

if (file.exists(here("output", "figures", "extended_data", "figure 1", "stan_output.rds"))) {

  # read saved model output
  samples <- readRDS(here("output", "figures", "extended_data", "figure 1", "stan_output.rds"))$samples
  x_pred <- readRDS(here("output", "figures", "extended_data", "figure 1", "stan_output.rds"))$x_pred
  
  data_read_from_file <- 1
  source(here("scripts", "model", "run_drift_diff.R"))
  
} else {
  
  data_read_from_file <- 0
  source(here("scripts", "model", "run_drift_diff.R"))
  
  # save model output
  stan_output <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(stan_output, file = here("output", "figures", "extended_data", "figure 1", "stan_output.rds"))
}

# create file to save output
path <- here("output", "figures", "extended_data", "figure 1")
dir.create(path)

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
if (length(neg_roots) > 1) {
  source(here("scripts", "model", "tipping_point.R"))
}

figure1 <- ggarrange(lake_men_plot + geom_hline(yintercept = most_prob_tp, linetype = "dashed", colour = "orange")  + theme(text = element_text(size = 7)),
                     NULL,
                     ts_plot_ordered + coord_flip(xlim = c(min(bga_data$X.stdlevel), max(bga_data$X.stdlevel))) +
                       theme(axis.ticks.x = element_blank(),
                             axis.title.y = element_blank()) +
                       geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") +
                       scale_y_discrete(labels = rep("", n_subjects)) + theme(text = element_text(size = 7)), 
                     ncol = 3, nrow = 1, widths = c(2, 0.2, 1))

drift_plot <- drift_plot + geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") + 
  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-0.05, 0.05)) + 
  theme(text = element_text(size = 7)) +
  geom_line(data = et_paper_drift, aes(x, y), size = 0.9, 
            colour = "red", linetype = "dashed")

drift_plot <- insert_xaxis_grob(drift_plot, CI_plot, grid::unit(.06, "null"), position = "top")

figure2 <- ggarrange(drift_plot,
                     stability_plot + theme(text = element_text(size = 7)),
                     ncol = 2, nrow = 1, align = "hv")

figure <- ggarrange(figure1, figure2, ncol = 1, nrow = 2, align = "hv", heights = c(1, 1))
figure

# save
ggsave(paste(path, "fig1.pdf", sep = "/"), figure, device = "pdf", height = 70, width = 88, units = "mm", dpi = 500)
