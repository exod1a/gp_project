
# code to make third figure in the publication
# roughly 22 mins to run first part

# Run cusp catastrophe process
# cusp parameters
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# set seed
seed <- 167 #231
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
# time step used for time series plot 
dt         <- 0.2
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

# create short time series plot
source(here("scripts", "model", "time_series_plot.R"))

ts_hist <- ggplot() + 
  geom_histogram(data = ts_data, aes(x=state), colour="black", fill="#f0d4ec", bins = 22) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "State", y = "Count") +
  theme_classic()

# create file to save output
path <- here("output", "figures", "main_text", "figure 3")
dir.create(path)

# process data
pred_inf <- 1
ground_truth <- 1
ts_data_full <- ts_data
source(here("scripts", "model", "prepare_data.R"))

if (file.exists(here("output", "figures", "main_text", "figure 3", "fig3.1_data.rds"))) {
  
  samples <- readRDS(here("output", "figures", "main_text", "figure 3", "fig3.1_data.rds"))$samples
  data_read_from_file <- 1
  source(here("scripts", "model", "run_drift_diff.R"))
} else {
  
  data_read_from_file <- 0
  source(here("scripts", "model", "run_drift_diff.R"))
  
  fig3_data <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(fig3_data, file = paste(path, "fig3.1_data.rds", sep = "/"))
}

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
source(here("scripts", "model", "tipping_point.R"))

drift_plot <- drift_plot + geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") + 
                           coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-2, 2)) + theme(text = element_text(size = 7))

drift_plot <- insert_xaxis_grob(drift_plot, CI_plot, grid::unit(.06, "null"), position = "top")

################################################################################################################################################################
# plot

#stability_df <- as.data.frame(table(multistability)/(ncol(drift_draws)-1))
#stability_df[5, 2] <- 0
#stability_df$multistability <- factor(stability_df$multistability, levels = c(0:4))
#stability_df[5, 1] <- 0

# plot
stability_plot <- ggplot(stability_df, 
                        aes(x = ms_subsetted, y = Freq)) +
  geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Multistability", y = "Probability") +
  theme_classic()

# plot tipping point region with drift and modality plot
figure1 <- plot_grid(ts_plot_ordered + geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") +
                       coord_cartesian(xlim = c(min(x_pred), max(x_pred))) + theme(text = element_text(size = 7)), 
                     ts_hist + coord_cartesian(xlim = c(min(x_pred), max(x_pred))) + theme(text = element_text(size = 7)),
                     drift_plot, 
                     stability_plot + theme(text = element_text(size = 7)),
                     diff_plot + geom_vline(xintercept = real_root, linetype = "dashed", colour = "orange") +
                       coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 1.5)) + theme(text = element_text(size = 7)), 
                     ncol = 2, nrow = 3, align = "hv")

# save
ggsave(paste(path, "fig3p1.pdf", sep = "/"), figure1, device = "pdf", height = 110, width = 90, units = "mm", dpi = 500)

################################################################################################################################################

# plot true positive heatmap for bistable case

true_pos_matrix <- matrix(NA, 4, 4)

t_c <- 49.04125

dt <- c(t_c, 0.1*t_c, 0.01*t_c, 0.001*t_c)
total_points <- c(25, 50, 100, 200)

#for (m in 1:length(dt)) {
#  for (n in 1:length(total_points)) {

for (m in 1:4) {
  for (n in 1:4) {
    
    true_positives <- read.table(here("output", "figures", "main_text", "figure 3", "heatmap_bistable", 
                                      paste(dt[m]/t_c, "t_c", "N", total_points[n], ".txt", sep = "")), header = FALSE, sep = "\n")
    
    true_pos_matrix[m, n] <- (sum(true_positives)) / (sum(true_positives) + sum(!true_positives)) * 100
  }
}

rownames(true_pos_matrix) <- c("t_c", "0.1t_c", "0.01t_c", "0.001t_c")
colnames(true_pos_matrix) <- c("25", "50", "100", "200")

to_plot <- melt(true_pos_matrix)
to_plot$Var1 <- factor(to_plot$Var1, ordered = T, levels = c("0.001t_c", "0.01t_c", "0.1t_c", "t_c"))
to_plot$Var2 <- factor(to_plot$Var2, ordered = T, levels = c("25", "50", "100", "200"))

bistable_heatmap <- ggplot(data = to_plot, aes(Var2, Var1, fill = value)) + 
  geom_tile() +
  #geom_text(aes(label=round(value, 0)), colour = "white") +
  scale_fill_gradientn(limits = c(0, 100), colours = c('midnightblue', '#FFFFFF', 'red')) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x = "Number of points", y = "Time step", title = "Bistable process, 100 simulations per cell", fill = "True positive rate (%)") +
  theme_classic() +
  theme(axis.line=element_blank(), text = element_text(size = 7))

ggsave(bistable_heatmap, filename = paste(path, "bistable_heatmap.pdf", sep = "/"), device = "pdf", dpi = 300,
       height = 45, width = 80, units = "mm")

################################################################################################################################################
# code to make the second part of the figure

#################################################################################################################

# Show that the model can distinguish between bimodality and bistability

# sample from a custom SDE density
# dens = vector of stationary density values 
# x = equally-spaced x-values associated with stationary density values
# N is the number of desired samples
# seed = seed
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

ground_truth <- 1
pred_inf <- 1
chosen_taxa <- 0

# Plotting
# Grid for x axis
lim    <- 1.4
x_grid <- seq(-lim, lim, length.out = 500)

roots <- x_grid[drift_roots(DRIFT(x_grid), x_grid)$root]

# custom drift plot
ground_truth_drift <- data.frame("x" = x_grid, "y" = DRIFT(x_grid))
drift_plot <- ggplot(ground_truth_drift) +
  geom_line(aes(x, y), size = 0.7) + labs(title = "", x = "State", y = "Drift") + 
  geom_hline(yintercept=0, size = 0.4) +
  geom_point(aes(roots, 0), size = 2.5) +
  coord_cartesian(xlim = c(min(x_grid), max(x_grid))) + 
  theme_classic()

# custom diffusion plot
ground_truth_diff <- data.frame("x" = x_grid, "y" = 0.5*DIFFUSION(x_grid))
diff_plot <- ggplot(ground_truth_diff) +
  geom_line(aes(x, y), size = 0.7) + labs(title = "", x = "State", y = "Diffusion") +
  coord_cartesian(xlim = c(min(x_grid), max(x_grid))) +
  theme_classic()

# custom density plot
ground_truth_stat_dens <- data.frame("x" = x_grid, 
                                     "y" = nonparametric_stationary_density(DRIFT(x_grid), 
                                                                            0.5*DIFFUSION(x_grid), 
                                                                            x_grid))
dens_plot <- ggplot(ground_truth_stat_dens) +
  labs(title = "", x = "State", y = "Density") +
  geom_ribbon(aes(x, ymin = 0, ymax = y), fill = "#f0d4ec", alpha = 1, colour = "black") +
  coord_cartesian(xlim = c(min(x_grid), max(x_grid))) + 
  theme_classic()

# sample short time series
# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01
# time step used for time series plot
dt         <- 1
num_points <- 2
# total number of points across all time series
total_points <- 100

if (dt * num_points >= 1) {
  end_time <- ceiling(dt * num_points)
} else {
  end_time <- round(dt * num_points, digits = 2)
}
time <- seq(from = 0, to = end_time, by = delta)


#if (dt * num_points >= 1) {
#  end_time <- ceiling(dt * (num_points - 1))
#} else {
#  end_time <- round(dt * (num_points - 1), digits = 3)
#}
#time <- seq(from = 0, to = end_time, by = delta)

# number of short time series 
n_subjects <- ceiling(total_points / num_points)

ts_data <- data.frame(matrix(nrow = 0, ncol = 3)) 

seed <- 15
set.seed(seed)

# simulate short time series from process
for (i in 1:n_subjects) {
  short_ts_data <- run_SDE(DRIFT, DIFFUSION, time, ground_truth_stat_dens) %>%
    cbind(time = time, subject = i) %>%
    as.data.frame() %>% 
    set_colnames(c("state", "time", "subject"))
  
  # subset data (otherwise it will take too long in Stan)
  ts_data <- rbind(ts_data, short_ts_data[seq(from = 1, to = length(time) - floor(dt/delta), by = floor(dt/delta)), ])
  
  #ts_data <- rbind(ts_data,
  #                 short_ts_data[seq(from = 1, to = length(time),
  #                                   by = min(floor(dt/delta), length(time)-1)), ])
}

# plot short time series
source(here("scripts", "model", "time_series_plot.R"))
ts_plot_ordered

ts_plot_ordered <- ts_plot_ordered + coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  labs(title = "") +
  scale_x_continuous(breaks = c(-1, 0, 1))
#ggsave(ts_plot_ordered, filename = here("data", "bimodal_ts.png"), dpi = 300,
#       height = 4, width = 6, units = "in")
# histogram of short time series data points
ts_hist <- ggplot() + 
  geom_histogram(data = ts_data, aes(x=state), colour="black", fill="#f0d4ec", binwidth = 0.09) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "State", y = "Count") +
  theme_classic()
ts_hist
#ggsave(ts_hist, filename = here("data", "bimodal_hist.png"), dpi = 300,
#       height = 4, width = 6, units = "in")
# create combined plot
SDE_plot <- ggarrange(ts_plot_ordered +
                        theme(axis.text.x = element_blank(),
                              axis.title.x = element_blank()) +
                        coord_cartesian(xlim = c(min(x_grid), max(x_grid))),
                      dens_plot + theme(axis.text.x = element_blank(),
                                        axis.title.x = element_blank()) +
                        coord_cartesian(xlim = c(min(x_grid), max(x_grid))),
                      drift_plot + theme(axis.text.x = element_blank(),
                                         axis.title.x = element_blank()) +
                        coord_cartesian(xlim = c(min(x_grid), max(x_grid))),
                      diff_plot + coord_cartesian(xlim = c(min(x_grid), max(x_grid))), 
                      ncol = 1, align = "hv")

# run model
source(here("scripts", "model", "prepare_data.R"))


if (file.exists(here("output", "figures", "main_text", "figure 3", "fig3.2_data.rds"))) {
  
  samples <- readRDS(here("output", "figures", "main_text", "figure 3", "fig3.2_data.rds"))$samples
  data_read_from_file <- 1
  source(here("scripts", "model", "run_drift_diff.R"))
} else {
  
  data_read_from_file <- 0
  source(here("scripts", "model", "run_drift_diff.R"))
  
  fig3p2_data <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(fig3p2_data, file = paste(path, "fig3.2_data.rds", sep = "/"))
}

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))

stability_df$ms_subsetted <- factor(stability_df$ms_subsetted, levels = c(1:5))
stability_df[4, ] <- c(4, 0)
stability_df[5, ] <- c(5, 0)

# plot
stability_plot <- ggplot(stability_df, 
                        aes(x = ms_subsetted, y = Freq)) +
  geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) + 
  labs(x = "Multistability", y = "Probability") +
  theme_classic()

figure2 <- ggarrange(ts_plot_ordered + coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
                       scale_x_continuous(breaks = c(-1, 0, 1)) + theme(text = element_text(size = 7)),
                     ts_hist + coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
                       scale_x_continuous(breaks = c(-1, 0, 1)) + theme(text = element_text(size = 7)),
                     drift_plot + coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-1, 0.5)) +
                       scale_x_continuous(breaks = c(-1, 0, 1)) + theme(text = element_text(size = 7)), 
                     stability_plot + theme(text = element_text(size = 7)),
                     diff_plot + coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(0, 0.6)) +
                       scale_x_continuous(breaks = c(-1, 0, 1)) + theme(text = element_text(size = 7)),
                     ncol = 2, nrow = 3, align = "hv")

# save
ggsave(paste(path, "fig3p2.pdf", sep = "/"), figure2, device = "pdf", height = 110, width = 90, units = "mm", dpi = 500)

##########################################################################################################################################################

# heatmap for bimodal case

t_c <- 13.53338

dt <- c(t_c, 0.1*t_c, 0.01*t_c, 0.001*t_c)
total_points <- c(25, 50, 100, 200)

# plot true positive heatmap for bistable case

true_pos_matrix <- matrix(NA, 4, 4)

#for (m in 1:length(dt)) {
#  for (n in 1:length(total_points)) {

for (m in 1:4) {
  for (n in 1:4) {
    true_positives <- read.table(here("output", "figures", "main_text", "figure 3", "bimodal_simulations", 
                                      paste("dt", dt[m]/t_c, "t_c", "N", total_points[n], sep = ""), "true_positives.txt"), header = FALSE, sep = "\n")
    
    true_pos_matrix[m, n] <- (sum(true_positives)) / (sum(true_positives) + sum(!true_positives)) * 100
  }
}

rownames(true_pos_matrix) <- c("t_c", "0.1t_c", "0.01t_c", "0.001t_c")
colnames(true_pos_matrix) <- c("25", "50", "100", "200")

# fix factors and levels
to_plot <- melt(true_pos_matrix)
to_plot$Var1 <- factor(to_plot$Var1, ordered = T, levels = c("0.001t_c", "0.01t_c", "0.1t_c", "t_c"))
to_plot$Var2 <- factor(to_plot$Var2, ordered = T, levels = c("25", "50", "100", "200"))

# plot matrix 
bimodal_heatmap <- ggplot(data = to_plot, aes(Var2, Var1, fill = value)) + 
  geom_tile() +
  #geom_text(aes(label=value), colour = "white") +
  scale_fill_gradientn(limits = c(0, 100), colours = c('midnightblue', '#FFFFFF', 'red')) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x = "Number of points", y = "Time step", title = "Bimodal process, 100 simulations per cell", fill = "True positive rate (%)") +
  theme_classic() +
  theme(axis.line=element_blank(), text = element_text(size = 7))

ggsave(bimodal_heatmap, filename = paste(path, "bimodal_heatmap.pdf", sep = "/"), dpi = 300, device = "pdf",
       height = 45, width = 80, units = "mm")

figure <- ggarrange(figure1, 
                    figure2, 
                    bistable_heatmap, 
                    bimodal_heatmap, 
                    nrow = 2, ncol = 2)

ggsave(figure, filename = paste(path, "fig3.pdf", sep = "/"), dpi = 300, device = "pdf",
       height = 180, width = 180, units = "mm")

