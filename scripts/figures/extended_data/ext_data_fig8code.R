
# top taxa included in the extended HITChip data set
top_taxa <- readLines(here("data", "ext_HITChip", "taxa.txt"))
top_taxa <- top_taxa[-which(top_taxa == "Prevotella_oralis_et_rel")]
# bistable taxa
bistable_taxa <- top_taxa[c(12, 32, 46, 62)]
# taxa found to be bistable but not all are in Lahti et al. 2014

# prepare extended HITChip data and save output
ts_data_full <- list()
for (m in 1:length(bistable_taxa)) {
  
  chosen_taxon <- bistable_taxa[m]
  
  ts_data_full[[m]] <- readRDS(here("data", "ext_HITChip", chosen_taxon, "ts_data.rds"))
  names(ts_data_full[[m]])[4] <- "state"
}

names(ts_data_full) <- bistable_taxa

plots <- lapply(1:4, function(i) {
  lapply(1:4, function(j) {
    
    ggplot() +
      geom_point(
        aes(
          x = ts_data_full[[i]]$state,
          y = ts_data_full[[j]]$state
        ),
        alpha = 0.5,
        size = 0.8
      ) +
      labs(
        x = bistable_taxa[i],
        y = bistable_taxa[j]
      ) +
      theme_classic() +
      theme(axis.title.x = element_text(size = 10),
            axis.title.y = element_text(size = 10))
  })
}) |> unlist(recursive = FALSE)

p <- wrap_plots(plots, ncol = 4)

# save
ggsave(paste(here("output", "figures", "extended_data", "figure 8"), "fig8.pdf", sep = "/"), 
       p, device = "pdf", height = 300, width = 300, units = "mm", dpi = 500)


