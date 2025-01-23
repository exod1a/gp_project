
# Create ext. data figure 4
# approximate run time: approx. 1 hour
# the output for the other tipping elements in Lahti et al. 2014

# top taxa included in the extended HITChip data set
top_taxa <- readLines(here("data", "ext_HITChip", "taxa.txt"))
# select tipping elements from Lahti et al. 2014
tipping_elements <- top_taxa[c(5, 46, 47, 28, 61, 62)]
# c(12, 32, 46, 47, 63) taxa found to be bistable but not all are in Lahti et al. 2014

plot_list <- list()

# prepare extended HITChip data and save output
for (m in 1:length(tipping_elements)) {

  chosen_taxa <- tipping_elements[m]
  # read data for taxa if model has already been previously run on it
  if (file.exists(here("output", "ext_HITChip", chosen_taxa, "stan_output.rds"))) {
    # read saved model output
    samples <- readRDS(here("output", "ext_HITChip", chosen_taxa, "stan_output.rds"))$samples
    x_pred <- readRDS(here("output", "ext_HITChip", chosen_taxa, "stan_output.rds"))$x_pred
    
    data_read_from_file <- 1
    pred_inf <- 1
    ground_truth <- 0
    source(here("scripts", "model", "run_drift_diff.R"))
    
    # read in full time series data without changes (but still with large time step points removed)
    ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
    names(ts_data_full)[4] <- "state"
    
    # read in time series data with changes e.g. dx, dt
    ts_data <- readRDS(here("data", "ext_HITChip", chosen_taxa, "taxon_data.rds"))

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
    
  } else { 
    source(here("scripts", "model", "run_model_on_taxa.R")) # run model
  }
  
  # find stability
  source(here("scripts", "model", "stability.R"))
  # if bistable
  if (length(pos_roots) > 0) {
    # get tipping region
    source(here("scripts", "model", "tipping_point.R"))
  
    # add tipping point to drift plot 
    #drift_plot <- drift_plot + geom_vline(xintercept = mean_tp, linetype = "dashed", colour = "orange") + 
    #  coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-0.1, 0.1)) + theme(text = element_text(size = 7))
  
    drift_plot <- insert_xaxis_grob(drift_plot, CI_plot, grid::unit(.06, "null"), position = "top")
  }

  if (length(pos_roots) == 0) {
    drift_plot <- drift_plot + coord_cartesian(xlim = c(min(x_pred), max(x_pred)), ylim = c(-0.1, 0.1))
  }
  
  # make sure all stability plots have the same range
  missing_probs <- c(1:4)[!(c(1:4) %in% stability_df$ms_subsetted)]
  temp_df <- as.data.frame(matrix(c(missing_probs, rep(0, 3)), nrow = length(missing_probs), ncol = 2))
  colnames(temp_df) <- c("ms_subsetted", "Freq")
  stability_df$ms_subsetted <- factor(stability_df$ms_subsetted, levels = c(1:4))
  stability_df <- rbind(stability_df, temp_df)
  
  # plot
  stability_plot <- ggplot(stability_df, 
                           aes(x = ms_subsetted, y = Freq)) +
    geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                       breaks = c(0, 0.25, 0.5, 0.75, 1),
                       limits = c(0, 1)) +
    labs(x = "Multistability", y = "Probability") +
    theme_classic()

  # create subfigure
  sub_figure <- ggarrange(ts_plot_ordered + labs(title = chosen_taxa) + coord_cartesian(xlim = c(min(x_pred), max(x_pred))),
                        drift_plot,
                        #diff_plot + coord_cartesian(xlim = c(min(x_pred), max(x_pred))),
                        stability_plot,
                        ncol = 1, nrow = 3, align = "hv")

  plot_list[[m]] <- sub_figure
}

# create and save figure
figure <- ggarrange(plot_list[[1]], 
                    plot_list[[2]], 
                    plot_list[[3]], 
                    plot_list[[4]], 
                    plot_list[[5]], 
                    plot_list[[6]],
                    ncol = 6, nrow = 1, align = "hv")

path <- here("output", "figures", "extended_data", "figure 4")
dir.create(path)

# save
ggsave(paste(path, "fig4.pdf", sep = "/"), figure, device = "pdf", height = 170, width = 300, units = "mm", dpi = 500)
