library(pheatmap)
library(grid)
library(dplyr)
library(limma)
library(ComplexHeatmap)
library(circlize)

data <- read.csv('/Path/To/Data/RNAseq/DEG.csv', header=TRUE)


## Filter for Genes with | Log2FC | > 1 
genes_logfc_1 <- data %>%
  filter((E_Hi.D_120_logFC > 1 & Z_Hi.D_125_logFC > 1) |
           (E_Hi.D_120_logFC < -1 & Z_Hi.D_125_logFC < -1)) %>%
  select(Gene, E_Hi.D_120_pval, E_Hi.D_120_logFC, Z_Hi.D_125_pval, Z_Hi.D_125_logFC)
## Find Genes with p value < 0.05
genes_pval <- data %>%
  filter(E_Hi.D_120_pval < 0.05 & Z_Hi.D_125_pval < 0.05) %>%
  select(Gene, E_Hi.D_120_pval, E_Hi.D_120_logFC, Z_Hi.D_125_pval, Z_Hi.D_125_logFC)


### Identify relevant gene sets for data. Example based on p value and Fold change cutoff
## genes_to_keep <- c(genes_logfc_1$Gene, genes_pval$Gene)

## Or, list genes by hand:
genes_to_keep <- c("Gene1", "Gene2", "Gene3")
cols_to_keep <- c("Condition1", "Condition2", "Condition3")


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


### Bubble Plot
# Load necessary libraries
library(dplyr)
library(ggplot2)

# Read in primary and secondary pathway enrichment data
main_data <- read.csv('/Path/To/Main_Data.csv', header = TRUE)
# annotation_data <- read.csv('/Path/To/Annotation_Data.csv', header = TRUE)

# Merge datasets based on shared pathway names
# merged_data <- merge(main_data, annotation_data, by = 'pathway_name', all.x = TRUE)

# Filter annotation data to keep only relevant pathways for plotting
filtered_annotation_data <- main_data[main_data$plot != "NO", ]
data <- filtered_annotation_data

# Sort data based on enrichment score (ascending or descending)
sort_type = 'desc'
if(sort_type == 'desc'){
  sorted_data <- data %>%
    arrange(desc(enrichment_score))
  
}
if(sort_type == 'asc'){
  sorted_data <- data %>%
    arrange(enrichment_score)
}


# Select top N unique pathways (e.g., 52 pathways = 104 rows if duplicated per group)
top_n_pathways <- 52
selected_data <- data.frame()
unique_pathway_count <- 0

for (i in 1:nrow(sorted_data)) {
  current_pathway <- sorted_data$pathway_name[i]
  
  if (!current_pathway %in% selected_data$pathway_name) {
    pathway_rows <- sorted_data %>% filter(pathway_name == current_pathway)
    selected_data <- rbind(selected_data, pathway_rows)
    unique_pathway_count <- unique_pathway_count + 1
  }
  
  if (unique_pathway_count == top_n_pathways) break
}

# Filter for directionality (if needed)
upregulated_only <- filtered_annotation_data[filtered_annotation_data$same_direction == "Both upregulated", ]
downregulated_only <- filtered_annotation_data[filtered_annotation_data$same_direction == "Both downregulated", ]

# Define which pathways to visualize
# pathways_to_plot <- unique(downregulated_only$pathway_name)
# Optionally add pathways to a character vector manually:
pathways_to_plot <- c("histone H3 acetylation", "histone acetylation", 
                      "histone acetyltransferase complex", 
                      "Antigen processing and presentation", 
                      "peptide antigen binding", 
                      "MHC protein complex")

# Load full dataset containing enrichment and p-values across all comparisons
bubble_data <- read.csv('/Path/To/Complete_Enrichment_Data.csv', header = TRUE)

# Subset for selected pathways
plot_data <- bubble_data %>% filter(pathway_name %in% pathways_to_plot)

# Calculate -log10(p-value) for bubble size
plot_data <- plot_data %>%
  mutate(log_pval = -log10(pval))

# Clean up pathway labels for display
plot_data$pathway_name <- gsub("_", " ", plot_data$pathway_name)
clean_pathway_levels <- gsub("_", " ", pathways_to_plot)
plot_data$pathway_name <- factor(plot_data$pathway_name, levels = clean_pathway_levels)

# Set factor levels for the groups (X-axis order)
plot_data$group <- factor(plot_data$group, 
                          levels = c("E_Hi-D_120", "Z_Hi-D_125", "V_Hi-D_133", "P_Hi-D_125"))

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
png("/Path/To/Save/BubblePlot.png", width = 5, height = 5, units = "in", res = 250)

# Create the bubble plot
ggplot(plot_data, aes(x = group, y = pathway_name, size = log_pval, color = enrichment_score)) +
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








