

# after running all the taxa on CSC
# top taxa included in the extended HITChip data set
delta_elpd <- as.numeric(readLines(here("output", "figures", "extended_data", "figure 2", "elpd.txt")))
length(delta_elpd)

# top taxa included in the extended HITChip data set
chosen_taxa <- readLines(here("data", "ext_HITChip", "taxa.txt"))

# remove Prevotella oralis. No longer want to include in analysis
to_remove <- which(chosen_taxa == "Prevotella_oralis_et_rel")

delta_elpd <- delta_elpd[-to_remove]

p <- ggplot() +
  geom_point(aes(1:62, delta_elpd)) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(x = "Taxon", y = "ELPD difference", title = "Model Comparison: Flexible Diffusion vs. Constant Diffusion",
       subtitle = "Both models applied to the 62 most prevalent HITChip taxa") + 
  theme_minimal()

ggsave(paste(here("output", "figures", "extended_data", "figure 2"), "fig2.pdf", sep = "/"), 
       p, device = "pdf", height = 100, width = 200, units = "mm", dpi = 500)













# top taxa included in the extended HITChip data set
chosen_taxa <- readLines(here("data", "ext_HITChip", "taxa.txt"))
  
k <- 1
chosen_taxon <- chosen_taxa[k]
  
# read in time series data with changes e.g. dx, dt
ts_data <- readRDS(here("data", "ext_HITChip", chosen_taxon, "taxon_data.rds"))
  
# rename columns
names(ts_data)[c(5, 6)] <- c("state", "dx")

t_c <- calculate_tc(ts_data)

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

# Stan data
n_chains <- 4
n_iter   <- 2000
  
# prepare data for stan
stan_data <- list(N_real = nrow(ts_data),
                    N_pred = length(x_pred),
                    x_real = ts_data$state,
                    x_pred = x_pred,
                    dx = ts_data$dx,
                    dt = ts_data$dt,
                    c = c,
                    pred_inf = pred_inf)
  
# compile model
EM_model <- stan_model(file = here("scripts", "model", "EM_model.stan"))
  
n_cores <- detectCores()
# set number of cores to be the same as the number of chains or the maximum number of cores is below n_chains
n_cores <- min(n_chains, n_cores)
samples <- sampling(EM_model, stan_data, chains = n_chains, iter = n_iter, cores = n_cores/2, seed = 234)

# compile model
EM_cd_model <- stan_model(file = here("scripts", "model", "EM_constant_diff.stan"))
samples_cd <- sampling(EM_cd_model, stan_data, chains = n_chains, iter = n_iter, cores = n_cores/2, seed = 234)
  
library(loo)
log_lik    <- extract_log_lik(samples, merge_chains = TRUE)
log_lik_cd <- extract_log_lik(samples_cd, merge_chains = TRUE)
  
aggregate_log_lik <- function(log_lik, series_id) {
                              sapply(unique(series_id), function(s) {
                              rowSums(log_lik[, series_id == s, drop = FALSE])})}
  
log_lik_series    <- aggregate_log_lik(log_lik, ts_data$subject)
log_lik_cd_series <- aggregate_log_lik(log_lik_cd, ts_data$subject)
  
loo_gp <- loo(log_lik_series)
loo_cd <- loo(log_lik_cd_series)
  
elpd_gp <- loo_gp$estimates["elpd_loo", "Estimate"]
elpd_cd <- loo_cd$estimates["elpd_loo", "Estimate"]
  
delta_elpd <- elpd_gp - elpd_cd



