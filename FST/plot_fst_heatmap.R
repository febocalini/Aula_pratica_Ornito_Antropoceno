plot_fst_heatmap <- function(fst_matrix, title = " Weir & Cockerham's 1984 FST") {
  library(ggplot2)
  library(RColorBrewer)
  
  fst_matrix <- round(as.matrix(fst_matrix), 4)
  pops <- rownames(fst_matrix)
  
  # Transformar matriz em data frame longo
  df <- expand.grid(pop1 = pops, pop2 = pops)
  df$fst <- as.vector(fst_matrix)
  
  ggplot(df, aes(x = pop1, y = pop2, fill = fst)) +
    geom_tile(color = "black") +
    geom_text(
      aes(label = ifelse(is.na(fst), "", sprintf("%.3f", fst))),
      fontface = "bold",
      size = 4
    ) +
    scale_fill_gradientn(
      colours = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100),
      limits = c(0, 1),
      na.value = "white"
    ) +
    theme_minimal(base_size = 14) +
    labs(x = "", y = "", fill = "Fst", title = title) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}
