
#######TEST#############################################################################

# drift-based approach to multistability
multistability <- c()
num_roots <- list()

# get the draws that give stationary solutions. "NS" here means "non-stationary"
for (i in 1:(ncol(drift_draws)-1)) {

  # for poor resolution x_pred grids, some errors may occur because my root solver can't resolve the roots
  # I don't know how to fix this. It slightly inflates the number of draws with 0 modes
  val <- try(drift_roots(drift_draws[, i], x_pred))
  if (inherits(val, "try-error")) {
    num_roots[[i]] <- data.frame("root" = 0, "sign" = 0)
    multistability[i] <- 0
    next
  }
  num_roots[[i]] <- val
  
  # check if solution is stationary i.e. number of stable equilibria = number of tipping points + 1
  if (all(is.na(val))) { # 0 roots
    num_roots[[i]] <- data.frame("root" = 0, "sign" = 0)
    multistability[i] <- 0
  } else if (sum(num_roots[[i]]$sign == -1) == (sum(num_roots[[i]]$sign == 1) + 1)) {
    multistability[i] <- sum(num_roots[[i]]$sign == -1)
  } else {
    multistability[i] <- "NS"
  }
}

# remove NAs and renormalise since they are nonstationary solutions
#multistability <- na.omit(multistability)
# remove non-stationary posterior draws
ms_subsetted <- multistability

if (0 %in% multistability) {
  to_remove <- which(ms_subsetted == 0)
  ms_subsetted <- ms_subsetted[-to_remove]
  #num_roots <- num_roots[-to_remove]
}
if ("NS" %in% multistability) {
  to_remove <- which(ms_subsetted == "NS")
  ms_subsetted <- ms_subsetted[-to_remove]
  #num_roots <- num_roots[-to_remove]
}

stability_df <- as.data.frame(table(ms_subsetted)/(length(ms_subsetted)))

# plot
stability_plot <- ggplot(stability_df, 
                         aes(x = ms_subsetted, y = Freq)) +
  geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Multistability", y = "Probability") +
  theme_classic()

# most likely number of stable states
n_states <- as.numeric(names(which(table(ms_subsetted) == max(table(ms_subsetted)))))

#######TEST#############################################################################


# drift-based approach to multistability
#multistability <- c()
#num_roots <- list()
#for (i in 1:(ncol(drift_draws)-1)) {
  # for poor resolution x_pred grids, some errors may occur because my root solver can't resolve the roots
  # I don't know how to fix this. It slightly inflates the number of draws with 0 modes
#  val <- try(drift_roots(drift_draws[, i], x_pred))
#  if (inherits(val, "try-error")) {
#    num_roots[[i]] <- 0
#    multistability[i] <- 0
#    next
#  }
#  num_roots[[i]] <- val
#  
#  if (exists("chosen_taxa")) {
    # the second stable mode is at the end of the range for prevotella so some posterior draws don't find it
    # and will be reported as unistable even though they have a tipping point which will make it seem like 
    # prevotella is unistable when that's not what the output really gives
    # If this is applied to other taxa with just one stable mode, it tends to exaggerate how likely they are to
    # be bistable because it might find a tipping point at the end range
#    if (grepl("Prevotella", chosen_taxa)) {
#      multistability[i] <- sum(num_roots[[i]]$sign == 1) + 1
#    } else {
#      multistability[i] <- sum(num_roots[[i]]$sign == -1)
#    }
#  } else {
#    multistability[i] <- sum(num_roots[[i]]$sign == -1)
#  }
#}

# remove NAs and renormalise since they are nonstationary solutions
#multistability <- na.omit(multistability)
#if (0 %in% multistability) {
#  multistability <- multistability[-which(multistability == 0)]
#}
# change any NAs (drift draws without any roots) to zeros
#multistability[which(is.na(multistability))] <- 0

#stability_df <- as.data.frame(table(multistability)/(length(multistability)))

# plot
#stability_plot <- ggplot(stability_df, 
#                         aes(x = multistability, y = Freq)) +
#  geom_bar(stat = "identity", fill = "#8DA0CB", colour = "black") +
#  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
#  labs(x = "Multistability", y = "Probability") +
#  theme_classic()

# most likely number of stable states
#n_states <- as.numeric(names(which(table(multistability) == max(table(multistability)))))

