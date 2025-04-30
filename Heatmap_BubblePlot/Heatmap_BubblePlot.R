# Load necessary libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(grid)
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(circlize)



### Bubble Plot

# Read in primary and secondary pathway enrichment data
pathway_data <- read.csv('./BubblePlot.csv', header = TRUE)


get_top_n_pathways <- function(data, score_col, sort_type = 'desc', top_n_pathways = 10) {
  # Sort the data based on the dynamic score column
  
  ## score_col is the column with your enrichment score (color for bubble plot)
  
  sorted_data <- if (sort_type == 'desc') {
    data %>% arrange(desc(.data[[score_col]]))
  } else if (sort_type == 'asc') {
    data %>% arrange(.data[[score_col]])
  } else {
    stop("sort_type must be 'asc' or 'desc'")
  }
  
  # Initialize selected data
  selected_data <- data.frame()
  unique_pathway_count <- 0
  
  # Loop to collect top unique pathways
  for (i in 1:nrow(sorted_data)) {
    current_pathway <- sorted_data$pathway_name[i]
    
    if (!current_pathway %in% selected_data$pathway_name) {
      pathway_rows <- sorted_data %>% filter(pathway_name == current_pathway)
      selected_data <- rbind(selected_data, pathway_rows)
      unique_pathway_count <- unique_pathway_count + 1
    }
    
    if (unique_pathway_count == top_n_pathways) break
  }
  
  return(selected_data)
}

top_pathways <- get_top_n_pathways(data = pathway_data, 'enrichment', sort_type = "desc", top_n_pathways = 10)


# Calculate -log10(p-value) for bubble size
plot_data <- top_pathways %>%
  mutate(log_pval = -log10(p_value))

# Clean up pathway labels for display - optional
plot_data$pathway_name <- gsub("_", " ", plot_data$pathway_name)
clean_pathway_levels <- gsub("_", " ", unique(plot_data$pathway_name))
plot_data$pathway_name <- factor(plot_data$pathway_name, levels = clean_pathway_levels)

# Set factor levels for the groups (X-axis order)
plot_data$group <- factor(plot_data$group, 
                          levels = c("condition_1", "condition_2", "condition_3"))

# Define custom theme for plot
custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = 7, margin = margin(b = 5)),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.3, 'cm'),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, size = 7, face = "bold"),
      axis.text.y = element_text(size = 9, face = "bold"),
      axis.title.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Save plot to PNG
png("./BubblePlot.png", width = 5, height = 5, units = "in", res = 250)

# Create the bubble plot
ggplot(plot_data, aes(x = group, y = pathway_name, size = log_pval, color = enrichment)) +
  geom_point() +
  scale_color_gradient2(midpoint = 0, low = "darkblue", mid = "grey", high = "#FF6347") +
  labs(
    title = "",
    y = "Pathways",
    color = "Enrichment Score",
    size = "p-value"
  ) +
  scale_size_continuous(
    range = c(0, 5),
    breaks = -log10(c(0.05, 0.01, 0.0001, 0.000001, 0.00000001, 0.0000000001)),
    labels = c(0.05, 0.01, 0.0001, 0.000001, 0.00000001, 0.0000000001)
  ) +
  custom_theme()

dev.off()

























## Heatmap



### Identify relevant gene sets for data.
### Example based on column from pathway analysis (GSEA or L2P)
gene_list <- function(column){
  build_vector <- c()
  for(genes in column){
    build_vector <- c(build_vector, strsplit(genes, split = " ")[[1]])
  }
  return(unique(build_vector))
}

## To use this function:
## genes <- gene_list(Heatmap_data$genes)


## Or, list genes by hand:
genes_to_keep <- c("Gene1", "Gene2", "Gene3")
cols_to_keep <- c("Condition1", "Condition2", "Condition3")



DEG <- read.csv('./Heatmap.csv', header=TRUE, row.names = 1)
scaled_data <- scale(DEG)

# Define heatmap color palette
heatmap_colors <- colorRampPalette(c("steelblue", "white", "red"))(100)

# Define group labels for columns (e.g., different treatment groups)
group_labels <- rep(c("condition_1", "condition_2", "condition_3", "condition_4"), each = 5)  # Adjust for your groups
annotation_col <- data.frame(Drug = factor(group_labels, levels = c("condition_1", "condition_2", "condition_3", "condition_4")))
drug_palette <- c("condition_1" = "#8E7CBB", "condition_2" = "#A3D9A5", "condition_3" = "#F4CDA5", "condition_4" = "#F8B3C4") 
annotation_colors <- list(Drug = drug_palette[names(drug_palette) %in% levels(annotation_col$Drug)])


ha <- HeatmapAnnotation(
  df = annotation_col,
  col = annotation_colors,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Condition",                   
    title_gp = gpar(fontsize = 8),         # legend header
    labels_gp = gpar(fontsize = 5.5)       # condition labels
  )
)
png("./Heatmap.png", width = 5, height = 5, units = "in", res = 300)
heatmap <- Heatmap(
  scaled_data, 
  col = heatmap_colors, 
  row_title = NULL, column_title = NULL, name = "Expression", 
  show_row_dend = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  top_annotation = ha,  # Column annotations
  heatmap_legend_param = list(title = "Expression Level", title_gp = gpar(fontsize = 8),  # annotation legend header 
                              labels_gp = gpar(fontsize = 5.5)  # legend labels
                              ),
  row_title_gp = gpar(fontsize = 0),     # Adjust row title font size
  column_names_gp = gpar(fontsize = 7.5),   # Adjust x-axis label font size
  row_names_gp = gpar(fontsize = 7.5)      # Adjust y-axis label font size
  
)

# Draw heatmap
draw(heatmap)
dev.off()


