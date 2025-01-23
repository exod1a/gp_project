
# code to make fourth figure in the publication
# prepare extended HITChip data and save output
chosen_taxa <- "Prevotella_melaninogenica_et_rel"
chosen_taxon <- "Prevotella_melaninogenica_et_rel"
#source(here("scripts", "HITChip", "process_HC_data.R"))

# read in prevotella data with changes e.g. dx, dt
ts_data  <- readRDS(here("data", "ext_HITChip", chosen_taxa, "taxon_data.rds"))

# rename columns
names(ts_data)[c(5, 6)] <- c("state", "dx")

# re-order
ts_data <- ts_data[order(ts_data$state), ]
ts_data$drift <- (ts_data$dx) / (ts_data$dt)
ts_data$diff <- (ts_data$dx^2) / (2*ts_data$dt)

# read in full time series data without changes
ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
names(ts_data_full)[4] <- "state"

to_keep <- sort(unique(c(which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample1)),
                         which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample2)))))
ts_data_full <- ts_data_full[to_keep, ]

# plot time series
ts_data_full$mean_val <- NA
for (i in unique(ts_data_full$subject)) {
  subject_vals <- which(ts_data_full$subject == i)
  ts_data_full$mean_val[subject_vals] <- (min(ts_data_full$state[subject_vals]) + max(ts_data_full$state[subject_vals]))/2
}

ts_data_full$mean_val <- factor(ts_data_full$mean_val)

# only take small dt's
ts_data <- subset(ts_data, dt < 45)

c <- min(ts_data$state) + (max(ts_data$state) - min(ts_data$state)) / 2  # centre of the data

# grid for predictions
x_pred <- seq(min(ts_data$state)-0.3, max(ts_data$state)+0.3, length.out = 100)

ground_truth <- 0
pred_inf <- 1

stan_data <- list(N_real = nrow(ts_data),
                  N_pred = length(x_pred),
                  x_real = ts_data$state,
                  x_pred = x_pred,
                  dx = ts_data$dx,
                  dt = ts_data$dt,
                  c = c,
                  pred_inf = pred_inf)

ts_data$subject <- unlist(ts_data$subject)

samples <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$samples
x_pred <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$x_pred
data_read_from_file <- 1
source(here("scripts", "model", "run_drift_diff.R"))

# create file to save output
path <- here("output", "figures", "main_text", "figure 4")
dir.create(path)
#fig4_data_prevotella <- list("samples" = samples, "x_pred" = x_pred)
#saveRDS(fig4_data_prevotella, file = paste(path, "fig4_data_prevotella.rds", sep = "/"))

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
if (length(pos_roots) > 0) {
  source(here("scripts", "model", "tipping_point.R"))
}

ts_plot_ordered <- ts_data_full %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "royalblue4") +
  geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.4) +
  labs(x = "System state", y = "Subject") +
  #geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") +
  theme_classic() + 
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())


# save prevotella versions of plots
ts_plot_prevotella <- ts_plot_ordered

# histogram of short time series data points
ts_hist_prevotella <- ggplot() + 
  geom_histogram(data = ts_data_full, aes(x=state), colour="black", fill="#f0d4ec", binwidth = 0.3) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "System state", y = "Count") +
  #geom_vline(xintercept = most_prob_tp, linetype = "dashed", colour = "orange") +
  coord_cartesian(xlim = c(min(x_pred), max(x_pred))) +
  theme_classic()
ts_hist_prevotella

stability_plot_prevotella <- stability_plot
# plot
#modality_df$x <- rep(1, 5)
#modality_plot_prevotella <- ggplot(modality_df, 
#                                   aes(x = x, y = Freq, fill = multimodality)) +
#                            geom_bar(position = "stack", stat = "identity", colour = "black") +
#                            scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
#                            labs(x = "Multistability", y = "Probability") +
#                            theme_classic()

# repeat analysis from above but with Dialister as the chosen taxon 
# prepare extended HITChip data and save output
chosen_taxa <- "Dialister"
chosen_taxon <- "Dialister"
#source(here("scripts", "HITChip", "process_HC_data.R"))

# read in prevotella data with changes e.g. dx, dt
ts_data <- readRDS(here("data", "ext_HITChip", chosen_taxa, "taxon_data.rds"))

# rename columns
names(ts_data)[c(5, 6)] <- c("state", "dx")

# re-order
ts_data <- ts_data[order(ts_data$state), ]
ts_data$drift <- (ts_data$dx) / (ts_data$dt)
ts_data$diff <- (ts_data$dx^2) / (2*ts_data$dt)

# read in full time series data without changes
ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
names(ts_data_full)[4] <- "state"

# subset for small time step points
to_keep <- sort(unique(c(which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample1)),
                         which(ts_data_full$sample %in% unlist(subset(ts_data, dt < 45)$sample2)))))
ts_data_full <- ts_data_full[to_keep, ]

# plot time series
ts_data_full$mean_val <- NA
for (i in unique(ts_data_full$subject)) {
  subject_vals <- which(ts_data_full$subject == i)
  ts_data_full$mean_val[subject_vals] <- (min(ts_data_full$state[subject_vals]) + max(ts_data_full$state[subject_vals]))/2
}

ts_data_full$mean_val <- factor(ts_data_full$mean_val)

ts_plot_ordered <- ts_data_full %>%
  ggplot() +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "royalblue4") +
  geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.4) +
  labs(x = "System State", y = "Subject") +
  theme_classic() + 
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

# only take small dt's
ts_data <- subset(ts_data, dt < 45)

c <- min(ts_data$state) + (max(ts_data$state) - min(ts_data$state)) / 2  # centre of the data

# grid for predictions
x_pred <- seq(min(ts_data$state)-0.3, max(ts_data$state)+0.3, length.out = 100)

stan_data <- list(N_real = nrow(ts_data),
                  N_pred = length(x_pred),
                  x_real = ts_data$state,
                  x_pred = x_pred,
                  dx = ts_data$dx,
                  dt = ts_data$dt,
                  c = c)

ground_truth <- 0
pred_inf <- 1

ts_data$subject <- unlist(ts_data$subject)
plot(ts_data$state, ts_data$drift)

set.seed(234)

samples <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$samples
x_pred <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$x_pred
data_read_from_file <- 1
source(here("scripts", "model", "run_drift_diff.R"))

# create file to save output
path <- here("output", "figures", "main_text", "figure 4")
dir.create(path)
#saveRDS(samples, file = paste(path, "fig4_data_dialister.rds", sep = "/"))

source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))

# save dialister versions of plots

# histogram of short time series data points
ts_hist_dialister <- ggplot() + 
  geom_histogram(data = ts_data_full, aes(x=state), colour="black", fill="#f0d4ec", binwidth = 0.3) + 
  #geom_line(data = ground_truth_stat_dens, aes(x, y), linetype = "dashed", size = 1) +
  labs(x = "System state", y = "Count") +
  theme_classic()
ts_hist_dialister

stability_df$ms_subsetted <- factor(stability_df$ms_subsetted, levels = c(1:4))
stability_df[4, ] <- c(4, 0)

stability_plot_dialister <- ggplot(stability_df, aes(x = ms_subsetted, y = Freq)) +
                           geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
                           scale_y_continuous(labels = scales::percent_format(accuracy = 1), 
                                              breaks = c(0, 0.2, 0.4, 0.6, 0.8), limits = c(0, 0.8)) +
                           labs(x = "Multistability", y = "Probability") +
                           theme_classic()

ts_plot_prevotella <- ts_plot_prevotella + labs(title = "Prevotella melaninogenica et rel.")
ts_plot_prevotella <- insert_xaxis_grob(ts_plot_prevotella, CI_plot, grid::unit(.06, "null"), position = "top")


########################################################################################################################

# PCoA plot

set.seed(155)

library(rmarkdown)
library(microbiome)
library(knitr)
library(dplyr)
library(tidyr)
library(earlywarnings)
library(phyloseq)

path <- here("data", "ext_HITChip")

#processed_data_folder <- path
#input_data_folder <- path
#phyloseq.file <- paste(processed_data_folder, "Atlas10k.Rds", sep = "/")
#species_data <- paste(processed_data_folder, "atlas_species.Rds", sep = "/")
#metadata_file_complete <- paste(input_data_folder, "final_atlas_metadata_full_nonpublic.Rds", sep = "/")

#meta.variables <- c("sample", "subject", "twin_id", "project", "sample_type", "dna_extraction_method", 
#                    "sex", "nationality", "bmi", "age", "age_group", "time", "health_status", "health_info", 
#                    "antibio", "medication", "probiotics", "age_cohort", "bmi_class")

# Read the data
#A <- readRDS(phyloseq.file)
#atlas.metadata.full <- readRDS(metadata_file_complete)

# Accepted projects; fecal RBB samples with no antibio/medication
#A.timeseries <- pick_timeseries(A, atlas.metadata.full) 

# Adult fecal samples with RBB extraction with no reported health issues
#A.baseline   <- pick_baseline(A)     

# Species abundance table
#species.rpa <- readRDS(species_data)

# Probe abundance table
#probe_abundances <- readRDS(paste(processed_data_folder, "atlas_oligo.Rds", sep = "/"))

# Compositional version for extra checks
#Ac <- transform(A, "compositional")

# Read Phylym level data
#P <- readRDS(paste(processed_data_folder, "atlas_phylum.Rds", sep = "/"))

# Various sample sets
#source(here("scripts", "HITChip", "samplesets.R")) # Check if all is correct here

#global.vars.probe <- names(meta(A))[grep("probe", names(meta(A)))]
#global.vars.genus <- names(meta(A))[grep("genus", names(meta(A)))]
#global.vars.numeric <- c("age", "bmi")
#global.vars.host <- c("sex","nationality","health_status","age_cohort","bmi_class")
#global.vars.discrete <- c("sample_type","dna_extraction_method","sex","nationality","age_group","health_status","age_cohort","bmi_class")

# --------------------

# Data for CR 30.1.2021; modifications 25.6.2023
#pseq <- subset_samples(A.timeseries, age_group == "adult" & health_status %in% c("healthy", "noncompromised"))
#pseq.baseline <- subset_samples(A.baseline, age_group == "adult" & health_status %in% c("healthy", "noncompromised"))

# get rid of points with no time stamp
#pseq <- ps_drop_incomplete(pseq, vars = "time")
# order so that subjects aren't scattered  
#otu_table(pseq) <-  otu_table(pseq)[, order(sample_data(pseq)[[2]])]

# need to re-order time points by subject
#new_indices <- c(1:dim(sample_data(pseq))[1])
#for (i in unique(sample_data(pseq)$subject)) {
  
 # subject_vals <- which(sample_data(pseq)$subject == i)
  
  #for (j in subject_vals[-length(subject_vals)]) {
   # if (sample_data(pseq)$time[j+1] > sample_data(pseq)$time[j]) {
  #    next
  #  } else {
      # find the new order of the indices 
  #    new_indices[subject_vals] <- subject_vals[order(sample_data(pseq)$time[subject_vals])]
  #    break
  #  }
  #}
#}
#otu_table(pseq) <-  otu_table(pseq)[, new_indices]
#sample_data(pseq)[, c(2, 12)]

# 128 subjects
#pseq <- subset_samples(pseq, bmi_class %in% c("lean", "underweight", "overweight", NA))

#saveRDS(pseq, here("data", "ord_plot_pseq_data.rds"))

# ordination with highlighted dominant genera
#pseq <- ps_calc_dominant(pseq, rank = "Genus")

#sample_data(pseq)$dominant_Genus[sample_data(pseq)$dominant_Genus != "Prevotella_melaninogenica_et_rel"] <- "Other"

# All of the data processing had to be commented out above since we can't share the data.
# This is the final version needed to create the ordination plot.
pseq <- readRDS(here("data", "ord_plot_pseq_data.rds"))
pseq.ord <- ordinate(pseq, "PCoA", "bray")
#p1 <- plot_ordination(pseq, pseq.ord, type = "samples", color = "dominant_Genus") + geom_point(size = 2, alpha = 1) + 
#  scale_colour_brewer(type = "qual", palette = "Paired") + theme(text = element_text(size = 7)) + 
#  theme_classic()
#p1

chosen_taxa <- "Prevotella_melaninogenica_et_rel"
samples <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$samples
x_pred <- readRDS(here("data", "ext_HITChip", chosen_taxa, "stan_output.rds"))$x_pred
data_read_from_file <- 1
pred_inf <- 1
ground_truth <- 0
source(here("scripts", "model", "run_drift_diff.R"))
source(here("scripts", "model", "stat_dens.R"))
source(here("scripts", "model", "stability.R"))
if (length(pos_roots) > 0) {
  source(here("scripts", "model", "tipping_point.R"))
}


# read in full time series data without changes
ts_data_full <- readRDS(here("data", "ext_HITChip", chosen_taxa, "ts_data.rds"))
names(ts_data_full)[4] <- "state"

state_val <- c()
subjects <- unique(sample_data(pseq)$subject)

# used to colour the line segments: find if from one point to the next a transtion occurred or not
for (i in 1:length(subjects)) {
  
  # get indices that correspond to a single subject
  subject_indices <- which(sample_data(pseq)$subject == subjects[i])
  
  for (j in 1:(length(subject_indices)-1)) {
    
    if (ts_data_full$state[subject_indices[j+1]] < most_prob_tp & ts_data_full$state[subject_indices[j]] < most_prob_tp) {
      # no transition from left
      state_val <- c(state_val, "low")
      
    } else if (ts_data_full$state[subject_indices[j+1]] > most_prob_tp & ts_data_full$state[subject_indices[j]] < most_prob_tp) {
      # transition from left
      state_val <- c(state_val, "transition")
      
    } else if (ts_data_full$state[subject_indices[j+1]] > most_prob_tp & ts_data_full$state[subject_indices[j]] > most_prob_tp) {
      # no transition from right
      state_val <- c(state_val, "high")
      
      
    } else if (ts_data_full$state[subject_indices[j+1]] < most_prob_tp & ts_data_full$state[subject_indices[j]] > most_prob_tp) {
      # transition from right 
      state_val <- c(state_val, "transition")
      
    }
  }
}

xstart  <- c()
ystart  <- c()
xfinish <- c()
yfinish <- c()
# get start and end points to draw lines between points 
for (i in 1:length(subjects)) {
  
  # get indices that correspond to a single subject
  subject_indices <- which(sample_data(pseq)$subject == subjects[i])
  
  # find starting and ending points for arrows
  xstart  <- c(xstart, pseq.ord$vectors[subject_indices[-length(subject_indices)], 1])
  ystart  <- c(ystart, pseq.ord$vectors[subject_indices[-length(subject_indices)], 2])
  xfinish <- c(xfinish, pseq.ord$vectors[subject_indices[-1], 1])
  yfinish <- c(yfinish, pseq.ord$vectors[subject_indices[-1], 2])
}

# get colours for points ie find which state each point is in
state_val_point <- c()
for (i in 1:nrow(sample_data(pseq))) {
  if (ts_data_full$state[i] < most_prob_tp) {
    
    state_val_point[i] <- "low"
  } else {
    
    state_val_point[i] <- "high"
  }
}

line_seg_data <- data.frame("xs" = xstart, "xf" = xfinish, "ys" = ystart, "yf" = yfinish, "state" = state_val)

pcoa_plot <- ggplot() +
  geom_point(aes(x = pseq.ord$vectors[, 1], 
                 y = pseq.ord$vectors[, 2], colour = state_val_point), alpha = 0.8) +
  geom_segment(data = dplyr::filter(line_seg_data, state == "high"),  
               aes(x = xs, y = ys, xend = xf, yend = yf, colour = state), alpha = 0.8, size = 0.15) +
  geom_segment(data = dplyr::filter(line_seg_data, state == "low"),  
               aes(x = xs, y = ys, xend = xf, yend = yf, colour = state), alpha = 0.8, size = 0.15) +
  geom_segment(data = dplyr::filter(line_seg_data, state == "transition"),  
               aes(x = xs, y = ys, xend = xf, yend = yf, colour = state), alpha = 0.8) +
  scale_colour_manual(values = c("#1F78B4", "#A6CEE3", "orchid1")) +
  labs(x = "Axis 1 [35.9%]", y = "Axis 2 [11.2%]", title = "PCoA: Human gut microbiome. 128 Participants, 2-6 time points each") +
  theme_classic()

ggsave(here("output", "figures", "main_text", "figure 4", "pcoa.pdf"), 
       pcoa_plot, device = "pdf", height = 100, width = 140, units = "mm", dpi = 500)

##########################################################################################


# create and save figure
figure1 <- ggarrange(ts_plot_prevotella,
                    ts_plot_ordered + labs(title = "Dialister"),
                    ts_hist_prevotella, 
                    ts_hist_dialister,
                    stability_plot_prevotella,
                    stability_plot_dialister,
                    ncol = 2, nrow = 3, align = "hv")

figure <- ggarrange(figure1, pcoa_plot, nrow = 1)

# add title
#figure <- annotate_figure(figure, top = text_grob("Fig. 4: Model output for taxa previously reported as bistable", face = "bold", size = 14))
figure
# save
path <- here("output", "figures", "main_text", "figure 4")
ggsave(paste(path, "fig4.pdf", sep = "/"), figure, device = "pdf", height = 100, width = 200, units = "mm", dpi = 500)
