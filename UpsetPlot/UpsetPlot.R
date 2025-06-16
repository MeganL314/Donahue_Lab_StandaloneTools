#### Load libraries
library(ComplexHeatmap)
library(ComplexUpset)
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(circlize)
library(tibble)
library(tidyr)
library(readxl)

#### Load DEG

DEG <- read_excel('./Path/To/data.xlsx')
DEG <- as.data.frame(DEG)
head(DEG)


rownames(DEG) <- DEG$Gene

suffix = '_logFC'
fc_cols <- grep(suffix, colnames(DEG), value = TRUE)
DEG[,fc_cols] <- sapply(DEG[ , fc_cols], as.numeric)

head(DEG[,fc_cols]) ## See that Infinity is now 'Inf' which R should recognize as numeric!


#### Filter data for genes w/ ADJUSTED p value < 0.05
#### For each union of contrasts, which genes are in opposite directions??
get_opposite_direction <- function(df, pval_suffix = "_adjpval", pval_threshold = 0.05, fc_suffix = "_logFC", 
                                   fc_threshold = 0.1, csv = "genes_direction") {
  
  df$gene <- rownames(df)
  pval_cols <- grep(pval_suffix, colnames(df), value = TRUE)
  fc_cols <- grep(fc_suffix, colnames(df), value = TRUE)
  
  contrasts <- sub(pval_suffix, "", pval_cols)
  
  result_list <- list()
  
  for (i in 1:nrow(df)) {
    gene <- df$gene[i]
    sig_contrasts <- c()
    logfc_vals <- c()
    
    # Which contrasts are significant (p-value and FC)
    for (contrast in contrasts) {
      pval <- df[i, paste0(contrast, pval_suffix)]
      logfc <- df[i, paste0(contrast, fc_suffix)]
      
      if (!is.na(pval) && pval < pval_threshold &&
          !is.na(logfc) && abs(logfc) >= fc_threshold) {
        sig_contrasts <- c(sig_contrasts, contrast)
        logfc_vals <- c(logfc_vals, logfc)
      }
    }
    
    if (length(sig_contrasts) == 0) next
    
    nonsig_contrasts <- setdiff(contrasts, sig_contrasts)
    is_only_significant <- all(
      df[i, paste0(nonsig_contrasts, pval_suffix)] >= pval_threshold |
        is.na(df[i, paste0(nonsig_contrasts, pval_suffix)])
    )
    
    if (!is_only_significant) next
    
    signs <- sign(logfc_vals)
    if (length(unique(signs[signs != 0])) > 1) {
      combo <- paste(sort(sig_contrasts), collapse = "_AND_")
      out_row <- setNames(c(gene, logfc_vals), c("gene", sig_contrasts))
      
      if (!combo %in% names(result_list)) result_list[[combo]] <- list()
      result_list[[combo]][[length(result_list[[combo]]) + 1]] <- out_row
    }
  }
  
  df_final <- list()
  
  outfile <- paste0("./Path/To/Output/", csv, ".csv")
  all_contrast <- file(outfile, open = "wt")
  
  for (combo in names(result_list)) {
    rows <- result_list[[combo]]
    all_cols <- unique(unlist(lapply(rows, names)))
    all_cols <- all_cols[!is.na(all_cols) & all_cols != ""]
    all_cols <- c("gene", sort(setdiff(all_cols, "gene")))
    
    row_list <- lapply(rows, function(x) {
      x <- as.list(x)
      row <- setNames(vector("list", length(all_cols)), all_cols)
      row[names(x)] <- x
      row
    })
    
    df_out <- unnest(data.frame(do.call(rbind, row_list), stringsAsFactors = FALSE), cols = all_cols)
    writeLines(paste("Contrast:", combo), all_contrast)
    write.table(df_out, all_contrast, sep = ",", row.names = FALSE, col.names = TRUE, na = "")
    writeLines("", all_contrast)
  }
  
  close(all_contrast)
}

get_opposite_direction(DEG, fc_threshold = 0.1, csv = "genes_direction_FC_test")






get_upset_data <- function(df, pval_suffix = "_adjpval", pval_threshold = 0.05, fc_suffix = "_logFC", fc_threshold = .01, 
                           output_csv = "upset_summary.csv") {
  
  # Input columns
  pval_cols <- grep(pval_suffix, colnames(df), value = TRUE)
  fc_cols   <- grep(fc_suffix,   colnames(df), value = TRUE)
  contrasts <- sub(pval_suffix, "", pval_cols)
  
  # Binary matrix based on p-value
  binary_mat <- as.data.frame(sapply(df[pval_cols], function(col) col < pval_threshold))
  
  # Apply fold change threshold if specified
  if (!is.null(fc_threshold) && length(fc_cols) > 0) {
    for (fc_col in fc_cols) {
      contrast_name <- sub(fc_suffix, "", fc_col)
      pval_col <- paste0(contrast_name, pval_suffix)
      if (pval_col %in% colnames(binary_mat)) {
        binary_mat[[pval_col]] <- binary_mat[[pval_col]] & abs(df[[fc_col]]) > fc_threshold
      }
    }
  }
  
  # Set row names and strip suffixes
  rownames(binary_mat) <- rownames(df)
  colnames(binary_mat) <- gsub(paste0(pval_suffix, "$"), "", colnames(binary_mat))
  
  # Generate direction using final binary matrix - after FC filter (if needed)
  # Output vector initialize
  fc_direction <- character(nrow(df))
  
  # Loop over each row
  for (i in seq_len(nrow(df))) {
    active_contrasts <- colnames(binary_mat)[which(binary_mat[i, ] == TRUE)]
    
    if (length(active_contrasts) == 0) {
      fc_direction[i] <- "none"
      print("None are TRUE")
      next
    }
    
    fc_values <- numeric(length(active_contrasts))
    print(fc_values)
    for (j in seq_along(active_contrasts)) {
      con <- active_contrasts[j]
      fc_col <- paste0(con, fc_suffix)
      val <- df[i, fc_col]
      fc_values[j] <- if (!is.na(val)) val else NA
    }
    
    fc_values <- fc_values[!is.na(fc_values)]
    print(paste("Row", i, "- fc_values:", paste(fc_values, collapse = ", ")))
    print(paste("Row", i, "- length of fc_values:", length(fc_values)))
    
    if (length(fc_values) == 0) {
      fc_direction[i] <- "none"
    } else if (all(fc_values > 0)) {
      fc_direction[i] <- "up"
    } else if (all(fc_values < 0)) {
      fc_direction[i] <- "down"
    } else {
      fc_direction[i] <- "opposite"
    }
  }
  
  
  # Intersection combination labels
  combo_labels <- apply(binary_mat, 1, function(x) {
    current <- names(which(x))
    if (length(current) == 0) {
      return("None")
    } else {
      return(paste(current, collapse = "+"))
    }
  })
  
  gene_sets <- split(rownames(binary_mat), combo_labels)
  combo_table <- table(combo_labels)
  
  # Gene matrix for CSV output
  max_len <- max(lengths(gene_sets))
  gene_df <- data.frame(lapply(gene_sets, function(x) {
    length(x) <- max_len
    return(x)
  }), check.names = FALSE)
  
  gene_df[is.na(gene_df)] <- " "
  write.csv(gene_df, output_csv, row.names = FALSE)
  
  # Attach direction
  binary_mat$direction <- fc_direction
  
  return(binary_mat)
}

#### Get Upset data for 
upset_binary_data <- get_upset_data(DEG, fc_threshold = 1, output_csv = "./Path/To/Output/Upset_Summary_p05_FC1_DEFindMarkers_Volcano_ID_test.csv")


## remove 'opposite'
upset_binary_data <- upset_binary_data[(upset_binary_data$direction != 'opposite'),]
upset_binary_data[is.na(upset_binary_data)] <- FALSE

## If you want to remove all rows where everything == FALSE (the gene is not significant for any contrast)
## Remove rows where everything is false
binary_mat_filtered <- upset_binary_data[rowSums(upset_binary_data[1:5]) > 0, ]

## Convert binary TRUE/FAlSE to 0/1 int
binary_mat_filtered[ , 1:(ncol(binary_mat_filtered) - 1)] <- lapply(binary_mat_filtered[ , 1:(ncol(binary_mat_filtered) - 1)], as.integer)



intersection_counts <- binary_mat_filtered %>% group_by(across(everything())) %>%
  summarise(count = n(), .groups = "drop") %>% arrange(desc(count))

write.csv(intersection_counts, "./Path/To/Output/intersection_counts.csv", row.names=FALSE)



contrast_order <- c("Example_10nM_vs_isotype_10nM", "Example_vs_isotype_1nM", "Example_vs_media",
                    "Example_10nM_vs_isotype_10nM", "Example_1nM_vs_isotype_1nM")


### For help with formating plots - changing font size, order, etc.:
### https://krassowski.github.io/complex-upset/articles/Examples_R.html#adjusting-set-size



#Plot the UpSet plot
png("./Path/To/Output/upset_plot_p05_FC_1_stackedbar.png", width = 8.5, height = 5, units = "in", res = 300)

ht <- upset(
  binary_mat_filtered,
  intersect = rev(contrast_order),
  name = "Significant Contrasts",
  sort_intersections=FALSE,
  intersections = list('Example_10nM','Example_10nM_vs_isotype_10nM', 'Example_vs_media',
                       'Example_1nM_vs_isotype_1nM', 'Example_1nM_vs_isotype_1nM',
                       c('Example_10nM_vs_isotype_10nM', 'Example_10nM_vs_isotype_10nM')),
  base_annotations = list(
    'Intersection size' = intersection_size(counts = TRUE, aes(fill = direction),
                                            text=list(size=2)) + ## adjust text size
      scale_fill_manual(values = c(up = "#1b9e77", down = "#d95f02", opposite = "#7570b3")) +
      theme(axis.title.y = element_text(hjust = .5, vjust = -30)
            )
  ),
  set_sizes = (
    upset_set_size() + 
      geom_text(aes(label = after_stat(count)),
      stat = "count", hjust = -0.1, size = 3, color = 'lightgrey') +
      theme(axis.text.x = element_text(angle = 90))
  ),
  width_ratio = 0.3,
  sort_sets = FALSE
)

print(ht)

dev.off()




