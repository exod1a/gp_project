
# make table of taxa that were studied for reviewer comments

# top taxa included in the extended HITChip data set
top_taxa <- readLines(here("data", "ext_HITChip", "taxa.txt"))
top_taxa <- top_taxa[-which(top_taxa == "Prevotella_oralis_et_rel")]
# bistable taxa
bistable_taxa <- top_taxa[c(12, 32, 46, 62)]

stability <- rep("unistable", length(top_taxa))
stability[c(12, 32, 46, 62)] <- "bistable"
exit_times <- rep("not calculated", length(top_taxa))
exit_times[c(12, 32, 46, 62)] <- "calculated"

df <- data.frame("System" = c("CUSP Model Simulations", "Custom SDE Model Simulations", "Lake Mendota"), 
                 "Stability" = c("bistable", "unistable", "bistable"), 
                 "Exit Time" = c("calculated", "not calculated", "calculated"))

df <- rbind(df, data.frame("System" = top_taxa, "Stability" = stability, "Exit Time" = exit_times))

df$type <- c(rep("other", 3), rep("taxon", 62))

library(gt)
library(scales)

n <- nrow(df)
mid <- ceiling(n / 2)

df_left  <- df[1:mid, ]
df_right <- df[(mid + 1):n, ]

# pad right side if needed
if (nrow(df_right) < nrow(df_left)) {
  df_right <- rbind(df_right,
                    data.frame(System = "", Stability = "", Exit.Time = "", type = ""))
}

df_wide <- data.frame(
  System_1 = df_left$System,
  Stability_1 = df_left$Stability,
  Exit_1 = df_left$Exit.Time,
  type_1 = df_left$type,
  System_2 = df_right$System,
  Stability_2 = df_right$Stability,
  Exit_2 = df_right$Exit.Time,
  type_2 = df_right$type)

taxa_table <- df_wide %>%
  gt() %>%
  tab_header(
    title = "The systems studied & the 62 most prevalent taxa selected for analysis"
  ) %>%
  cols_label(
    System_1 = "System",
    Stability_1 = "Stability",
    Exit_1 = "Exit Time",
    System_2 = "System",
    Stability_2 = "Stability",
    Exit_2 = "Exit Time"
  ) %>%
  tab_spanner(
    label = " ",
    id = "left",
    columns = c(System_1, Stability_1, Exit_1)
  ) %>%
  tab_spanner(
    label = " ",
    id = "right",
    columns = c(System_2, Stability_2, Exit_2)
  )

taxa_table <- taxa_table %>%
  data_color(
    columns = c(Stability_1, Stability_2),
    colors = scales::col_factor(
      palette = c("bistable" = "#a6dba0", "unistable" = "#fafafa"),
      levels = c("bistable", "unistable")
    )
  ) %>%
  data_color(
    columns = c(Exit_1, Exit_2),
    colors = scales::col_factor(
      palette = c("calculated" = "#a6dba0", "not calculated" = "#fafafa"),
      levels = c("calculated", "not calculated")
    )
  )

taxa_table <- taxa_table %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  )

taxa_table <- taxa_table %>%
  tab_style(
    style = cell_text(weight = "bold", size = px(28)),
    locations = cells_title(groups = "title")
  )

taxa_table <- taxa_table %>%
  cols_hide(columns = c(type_1, type_2)) %>%
  tab_style(
    style = cell_fill(color = "lavenderblush2"),
    locations = cells_body(
      columns = System_1,
      rows = type_1 == "other"
    )
  ) 

taxa_table <- taxa_table %>%
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "black",
      weight = px(2)
    ),
    locations = cells_body(columns = Exit_1)
  ) %>%
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "black",
      weight = px(2)
    ),
    locations = cells_column_labels(columns = Exit_1)
  )

library(webshot2)

# create file to save output
path <- here("output", "figures", "extended_data", "figure 6")
dir.create(path)

# Save as PDF
gtsave(taxa_table, here("output", "figures", "extended_data", "figure 6", "fig6.png"))

