
# Run model on extended HITChip taxa, subsetting for dt < 45 (high quality time points)

# pass chosen_taxa list to this file before running
# make sure "taxon_data.rds" and "ts_data.rds" files already exist as well
# which come from running "model/HITChip/process_HC_data.R"

for (k in 1:length(chosen_taxa)) {
  chosen_taxon <- chosen_taxa[k]
  #source(here("scripts", "HITChip", "process_HC_data.R"))  
  
  # read in time series data with changes e.g. dx, dt
  ts_data <- readRDS(here("data", "ext_HITChip", chosen_taxon, "taxon_data.rds"))
  
  # rename columns
  names(ts_data)[c(5, 6)] <- c("state", "dx")
  
  # re-order
  ts_data <- ts_data[order(ts_data$state), ]
  ts_data$drift <- (ts_data$dx) / (ts_data$dt)
  ts_data$diff <- (ts_data$dx^2) / (2*ts_data$dt)
  
  # read in full time series data without changes (but still with large time step points removed)
  ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxon, "ts_data.rds"))
  names(ts_data_full)[4] <- "state"
  
  # subset for small time step points
  to_keep <- sort(unique(c(which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample1)),
                           which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample2)))))
  ts_data_full <- ts_data_full[to_keep, ]
  
  # plot time series
  ts_data_full$mean_val <- NA
  for (j in unique(ts_data_full$subject)) {
    subject_vals <- which(ts_data_full$subject == j)
    ts_data_full$mean_val[subject_vals] <- (min(ts_data_full$state[subject_vals]) + max(ts_data_full$state[subject_vals]))/2
  }
  
  # add variable for ordering short time series
  ts_data_full$mean_val <- factor(ts_data_full$mean_val)

  # plot short time series
  ts_plot_ordered <- ts_data_full %>%
    ggplot() +
    geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "royalblue4") +
    geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.4) +
    labs(x = "System State", y = "Subject") +
    theme_classic() + 
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  
  # only take small dt's
  ts_data <- subset(ts_data, dt < 45)
  
  # find centre of data
  c <- min(ts_data$state) + (max(ts_data$state) - min(ts_data$state)) / 2  # centre of the data
  
  # grid for predictions
  x_pred <- seq(min(ts_data$state) - 0.3, max(ts_data$state) + 0.3, length.out = 100)
  
  # to make sure "run_drift_diff.R" runs
  data_read_from_file <- 0
  # we want predictive inference in the Gaussian process
  pred_inf <- 1
  # there is no ground truth function
  ground_truth <- 0
  
  # prepare data for stan
  stan_data <- list(N_real = nrow(ts_data),
                    N_pred = length(x_pred),
                    x_real = ts_data$state,
                    x_pred = x_pred,
                    dx = ts_data$dx,
                    dt = ts_data$dt,
                    c = c,
                    pred_inf = pred_inf)
  
  set.seed(155)
  # Stan data
  n_chains <- 4
  n_iter   <- 2000
  source(here("scripts", "model", "run_drift_diff.R"))
  
  # screate path for output
  path <- here("output", "ext_HITChip", chosen_taxon)
  dir.create(path)
  # save model output
  stan_output <- list("samples" = samples, "x_pred" = x_pred)
  saveRDS(stan_output, file = paste(path, "stan_output.rds", sep = "/"))
}

