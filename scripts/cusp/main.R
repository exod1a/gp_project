
# See https://www.geogebra.org/m/erbtjvzt for an interactive demo to get an understanding of the role of the parameters.
#  The demo plots the potential landscape in which the particle travels and the corresponding stationary density

# Run cusp catastrophe process
# cusp parameters
r       <- 0.2
alpha   <- -0.5
beta    <- 6
lambda  <- 4
epsilon <- 1.2

# Time grid. Large resolution (by = 0.01) is needed to obtain convergent simulation
delta      <- 0.01  # time step
end_time   <- 500   # number of integer points
time       <- seq(from = 0, to = end_time - 1, by = delta)
# time step used for time series plot
dt         <- 1
num_points <- ceiling(length(time) / (dt/delta))

seed <- 112
set.seed(seed)

# create file to save output
path <- here("output", "cusp", paste("r", r, "a", alpha, "b", beta, "eps", epsilon, "seed", seed, 
                                                "time", end_time, "dt", dt, sep = "_"))
dir.create(path)

source(here("scripts", "cusp", "run_cusp.R"))
saveRDS(ts_data, file = here(path, "ts_data.rds"))  # save data
source(here("scripts", "cusp", "plots.R"))
ggsave(here(path, "cusp_plot.pdf"), cusp_plot, height = 12, width = 8, units = "in", dpi = 300)  # save output
