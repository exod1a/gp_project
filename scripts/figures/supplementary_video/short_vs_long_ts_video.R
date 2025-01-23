
# Make supplementary video comparing long and short time series

source(here("scripts", "figures", "main_text", "short_vs_long_ts.R"))

# try to make a video showing this

subsetted_time <- 250

lts <- df_lts[1:(2*subsetted_time+1), c(1, 4)]
lts$step <- rep(c(0:subsetted_time), each = 2)[-(subsetted_time*2+2)]
colnames(lts)[2] <- "state"

# long time series animation
lts_anim <- ggplot(lts) +
  geom_line(aes(x=time, y=state), colour = "royalblue4", size = 0.2) +
  geom_point(aes(x=time, y=state), colour = "royalblue4", size = NA) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "Time", y = "System state", title = "Long time series") +
  coord_flip() +
  ylim(c(0, 8)) +
  theme_classic() +
  transition_reveal(step)

# break into frames
lts_ani <- lts %>% 
  split(.$step) %>% 
  accumulate(~ bind_rows(.x, .y)) %>% 
  bind_rows(.id = "frame") %>% 
  mutate(frame = as.integer(frame))

# long time series histogram animation
lts_hist_anim <- ggplot(lts_ani) +
  geom_histogram(aes(x=state, y = ..density..), fill = "slateblue", colour = "black") +
  geom_line(data = ground_truth_stat_dens, aes(x, y), size = 1) + 
  theme_classic() +
  ylim(c(0, 1)) +
  xlim(c(0, 8)) +
  theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
  transition_manual(frame)

lts_gif <- animate(lts_anim, end_pause = 30, fps = 10, width = 600, height = 600, res = 300)
lts_hist_gif <- animate(lts_hist_anim, end_pause = 30, fps = 10, width = 600, height = 600, res = 300)

# combine animations
lts_combined <- image_append(c(lts_hist_gif[1], lts_gif[1]), stack = T)
for(i in 2:100){
  combined <- image_append(c(lts_hist_gif[i], lts_gif[i]), stack = T)
  lts_combined <- c(lts_combined, combined)
}

#lts_combined

#############################################

sts <- df_sts[1:(2*subsetted_time), c(1, 2)]
sts$subject <- rep(c(1:subsetted_time), each = 2)
colnames(sts)[2] <- "state"

sts$mean_val <- NA
for (i in unique(sts$subject)) {
  subject_vals <- which(sts$subject == i)
  sts$mean_val[subject_vals] <- (min(sts$state[subject_vals]) + max(sts$state[subject_vals]))/2
}

sts$mean_val <- factor(sts$mean_val)

# short time series animation
sts_anim <- ggplot(sts) +
  geom_line(aes(x=state, y=mean_val, group=as.character(mean_val)), colour = "royalblue4", size = 0.2) +
  #geom_point(aes(x = state, y = mean_val, group = as.character(mean_val)), size = 0.2) +
  labs(x = "System state", y = "Subject", title = "Short time series") +
  theme_classic() + 
  xlim(c(0, 8)) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  transition_reveal(subject)

# break into frames
sts_ani <- sts %>% 
  split(.$subject) %>% 
  accumulate(~ bind_rows(.x, .y)) %>% 
  bind_rows(.id = "frame") %>% 
  mutate(frame = as.integer(frame))

# short time series histogram animation 

sts_hist_anim <- ggplot(sts_ani) +
  geom_histogram(aes(x=state, y = ..density..), fill = "slateblue", colour = "black") +
  geom_line(data = ground_truth_stat_dens, aes(x, y), size = 1) + 
  theme_classic() +
  ylim(c(0, 1)) +
  xlim(c(0, 8)) +
  theme(axis.line.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
        axis.line.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank()) +
  transition_manual(frame)

sts_gif <- animate(sts_anim, end_pause = 30, fps = 10, width = 600, height = 600, res = 300)
sts_hist_gif <- animate(sts_hist_anim, end_pause = 30, fps = 10, width = 600, height = 600, res = 300)

# combine animations
sts_combined <- image_append(c(sts_hist_gif[1], sts_gif[1]), stack = T)
for(i in 2:100){
  combined <- image_append(c(sts_hist_gif[i], sts_gif[i]), stack = T)
  sts_combined <- c(sts_combined, combined)
}

#sts_combined

#############################################

# combine all  animations
ts_gifs <- image_append(c(sts_combined[1], lts_combined[1]))
for (i in 2:100){
  combined <- image_append(c(sts_combined[i], lts_combined[i]))
  ts_gifs <- c(ts_gifs, combined)
}

#ts_gifs

## add error plot

error_lts <- c()
error_sts <- c()

# calculate density errors over all 100 long time series
for (i in 1:subsetted_time) {
  
  # find where to cut off time series for a given amount of time 
  time_bound <- min(which(as.integer(lts[,1]) == i))
  
  # density calculation
  approx_dens_lts <- density(lts[1:time_bound, 2], 
                             n = nrow(ground_truth_stat_dens), 
                             from = data_range[1], 
                             to = data_range[length(data_range)])
  
  # add a small number to any zeroes in the distribution in order to calculate the KL divergence
  approx_dens_lts$y[which(approx_dens_lts$y == 0)] <- approx_dens_lts$y[which(approx_dens_lts$y == 0)] + 
    .Machine$double.eps
  
  # long time series error
  error_lts[[i]] <- KL_div(ground_truth_stat_dens$y, approx_dens_lts$y)
}

# calculate density errors over all 100 sets of short time series
# calculate density errors over all 100 long time series
for (i in 1:subsetted_time) {
  
  # find where to cut off time series for a given amount of time 
  time_bound <- grep(1, sts[,1])[i]
  
  # density calculation
  approx_dens_sts <- density(sts[1:time_bound, 2], 
                             n = nrow(ground_truth_stat_dens), 
                             from = data_range[1], 
                             to = data_range[length(data_range)])
  
  # add a small number to any zeroes in the distribution in order to calculate the KL divergence
  approx_dens_sts$y[which(approx_dens_sts$y == 0)] <- approx_dens_sts$y[which(approx_dens_sts$y == 0)] + 
    .Machine$double.eps
  
  # short time series error
  error_sts[[i]] <- KL_div(ground_truth_stat_dens$y, approx_dens_sts$y)
}

# normalise data
alpha_lts <- 1 / max(unlist(error_lts))
alpha_sts <- 1 / max(unlist(error_sts))

adj_lts_err <- alpha_lts * unlist(error_lts)
adj_sts_err <- alpha_sts * unlist(error_sts)
# make it increase towards 1 instead of decrease towards 0
adj_lts_err <- 1 - adj_lts_err
adj_sts_err <- 1 - adj_sts_err

# make dataframe for labelling plot
error_data <- data.frame("time" = rep(1:subsetted_time, 2), 
                         "ts" = c(adj_lts_err, adj_sts_err), 
                         "type" = c(rep("Long time series", subsetted_time), rep("Short time series", subsetted_time)))

error_data$type <- factor(error_data$type, levels = c("Short time series", "Long time series"))

error_anim <- ggplot(error_data) +
  geom_line(aes(x = time, y = ts, colour = type, size = type, linetype = type)) +
  scale_colour_manual(values = c("black", "gray")) +
  scale_size_manual(values = c(0.7, 0.9)) +
  scale_linetype_manual(values = c("solid", "solid")) +
  labs(x = "Time", y = "") +
  ggtitle("<span style='font-size: 9pt;'>Agreement between data histogram and true stationary density</font>") +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "bottom", plot.title = element_markdown()) +
  transition_reveal(time)

error_gif <- animate(error_anim, end_pause = 30, fps = 10, width = 1200, height = 800, res = 300)

# combine all  animations
whole_gif <- image_append(c(ts_gifs[1], error_gif[1]), stack = T)
for (i in 2:100){
  combined <- image_append(c(ts_gifs[i], error_gif[i]), stack = T)
  whole_gif <- c(whole_gif, combined)
}

anim_save(animation = whole_gif, filename = here("output", "videos", "ts_comparison_fast_high_res.gif"))

