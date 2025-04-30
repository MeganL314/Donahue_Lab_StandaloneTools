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
## genes <- gene_list(data$genes)




## Or, list genes by hand:
genes_to_keep <- c("Gene1", "Gene2", "Gene3")
cols_to_keep <- c("Condition1", "Condition2", "Condition3")





data <- read.csv('/Path/To/Data/RNAseq/DEG.csv', header=TRUE)


df <- read.csv('/Path/to/your/data/.csv', header=TRUE, row.names = 1)

# Preprocess data (apply any transformations, e.g., log2 FC, scaling)
numeric_columns <- df[sapply(df, is.numeric)]
min_value <- min(numeric_columns, na.rm = TRUE)
shift_value <- if (min_value < 0) abs(min_value) + 1 else 0
pseudocount <- 1e-50

pairs_list <- list(
  c("TU1_E_Hi", "TU1_D_120"),
  c("TU2_E_Hi", "TU2_D_120"),
  c("TU3_E_Hi", "TU3_D_120"),
  c("TU4_E_Hi", "TU4_D_120"),
  c("TU5_E_Hi", "TU5_D_120"),
  c("TU1_Z_Hi", "TU1_D_125"),
  c("TU2_Z_Hi", "TU2_D_125"),
  c("TU3_Z_Hi", "TU3_D_125"),
  c("TU4_Z_Hi", "TU4_D_125"),
  c("TU5_Z_Hi", "TU5_D_125"))

# Assuming you are calculating Log2 Fold Change or similar metrics
for (col_pair in pairs_list) {
    col1 <- col_pair[1]
    col2 <- col_pair[2]
  # New column name based on pairs data
    col1_parts <- strsplit(col1, "_")[[1]]
    col2_parts <- strsplit(col2, "_")[[1]]
    new_column_name <- paste(col1_parts[1], col1_parts[2], "vs", col2_parts[2], col2_parts[3], sep = "_")
    # Compute the log2 fold change
    df[[new_column_name]] <- log2((df[[col1]] + pseudocount + shift_value) / (df[[col2]] + pseudocount + shift_value))
}

df_subset <- df[genes_to_keep, cols_to_keep, drop = FALSE]
scaled_data <- scale(df_subset)

# Define heatmap color palette
heatmap_colors <- colorRampPalette(c("steelblue", "white", "red"))(100)

# Define group labels for columns (e.g., different treatment groups)
group_labels <- rep(c("Ent", "Zaba"), each = 5)  # Adjust for your groups
heatmap_colors <- colorRampPalette(c("steelblue", "white", "red"))(100)
annotation_col <- data.frame(Drug = factor(group_labels, levels = c("Ent", "Zaba", "Vor", "Pano")))
drug_palette <- c("Ent" = "#8E7CBB", "Zaba" = "#A3D9A5", "Vor" = "#F4CDA5", "Pano" = "#F8B3C4") 
annotation_colors <- list(Drug = drug_palette[names(drug_palette) %in% levels(annotation_col$Drug)])

ha <- HeatmapAnnotation(df = annotation_col, col = annotation_colors)
png("/path/to/save/heatmap.png", width = 8, height = 6, units = "in", res = 300)

heatmap <- Heatmap(
  scaled_data, 
  col = heatmap_colors, 
  row_title = NULL, column_title = NULL, name = "Expression", 
  show_row_dend = FALSE, cluster_columns = FALSE, cluster_rows = FALSE,
  show_row_names = TRUE, show_column_names = TRUE, 
  top_annotation = ha,  # Column annotations
  heatmap_legend_param = list(title = "Expression Level")
)

# Draw heatmap
draw(heatmap)

# Close the PNG device
dev.off()


