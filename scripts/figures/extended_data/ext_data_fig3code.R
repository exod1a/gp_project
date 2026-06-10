
# plot true positive heatmap for bistable case
# plot the distance to the tipping point
true_pos_matrix <- matrix(NA, 4, 4)
tp_dist <- matrix(NA, 4, 4)

t_c <- 49.04125

dt <- c(t_c, 0.1*t_c, 0.01*t_c, 0.001*t_c)
total_points <- c(26, 50, 100, 200)

for (m in 1:length(dt)) {
  for (n in 1:length(total_points)) {
    
    heatmap_data <- read.table(here("output", "figures", "extended_data", "figure 3", 
                                    paste(dt[m]/t_c, "t_c", "N", total_points[n], ".txt", sep = "")), header = T)
    
    true_pos_matrix[m, n] <- (sum(heatmap_data$bistable)) / (sum(heatmap_data$bistable) + sum(!heatmap_data$bistable)) * 100
    
    tp_dist[m, n] <- mean(heatmap_data$tp_dist[!is.na(heatmap_data$tp_dist)])
  }
}

rownames(true_pos_matrix) <- c("t_c", "0.1t_c", "0.01t_c", "0.001t_c")
colnames(true_pos_matrix) <- c("13", "25", "50", "100")

rownames(tp_dist) <- c("t_c", "0.1t_c", "0.01t_c", "0.001t_c")
colnames(tp_dist) <- c("13", "25", "50", "100")

to_plot <- melt(true_pos_matrix)
to_plot$Var1 <- factor(to_plot$Var1, ordered = T, levels = c("0.001t_c", "0.01t_c", "0.1t_c", "t_c"))
to_plot$Var2 <- factor(to_plot$Var2, ordered = T, levels = c("13", "25", "50", "100"))

tp_plot <- melt(tp_dist)
tp_plot$Var1 <- factor(tp_plot$Var1, ordered = T, levels = c("0.001t_c", "0.01t_c", "0.1t_c", "t_c"))
tp_plot$Var2 <- factor(tp_plot$Var2, ordered = T, levels = c("13", "25", "50", "100"))

tp_error_heatmap <- ggplot(data = tp_plot, aes(Var2, Var1, fill = 1/value)) + 
  geom_tile() +
  geom_text(data = to_plot, aes(x = Var2, y = Var1, label = round(value, 0)), colour = "black") +
  scale_fill_gradientn(limits = c(0.8, 1.55), colours = c('midnightblue', '#FFFFFF', 'red')) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(x = "Number of time series", y = "Time step", title = "Tipping point esimation accuracy", 
       subtitle = "Text: true positive rate (%)", fill = "Inverse distance to tipping pt \n conditional on detection") +
  theme_classic() +
  theme(axis.line=element_blank(), text = element_text(size = 7))


#ggplot(tp_plot, aes(x = Var2, y = value, group = Var1, colour = Var1)) +
#  geom_line(size = 1) +
#  geom_point() +
#  labs(x = "Number of time series", y = "Distance to tipping point\n conditional on detection", title = "Tipping point esimation accuracy", 
#       colour = "Time step") + 
#  theme_minimal()


ggsave(tp_error_heatmap, filename = here("output", "figures", "extended_data", "figure 3", "heatmap.pdf"), device = "pdf", dpi = 300,
       height = 80, width = 110, units = "mm")
