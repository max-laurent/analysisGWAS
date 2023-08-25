LD_heatmap_plot <- function(data, x = "index_1", y = "index_2", locus1 = "locus1", locus2 = "locus2", LDcol = "r2k", limits = 100){


  data <- data %>% filter(., .data[[x]] < 100 & .data[[y]] <100 & .data[[locus1]] != .data[[locus2]])
  # Create a list containing the information for the mapping of the data
  mapping <- aes(x = .data[[x]], y= .data[[y]], fill = .data[[LDcol]])


  ggplot(data, mapping = mapping) +
    geom_tile() +
    scale_fill_gradient(low = "yellow", high = "red") +
    theme_bw()
}
