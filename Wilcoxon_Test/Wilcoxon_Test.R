library(tidyverse)
library(dplyr)
# import raw data
# to create this, use the import feature in R studio

library(readxl)
rawdata <- read_excel("./Wilcoxon_Test/.xlsx", sheet = "Fake_Paired_Analyte_Data")


custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=13, face = "bold", margin = margin(t = 0, r = 0, b = 8, l = 0)),
      legend.position = "right",
      legend.title = element_blank(),
      axis.text.y = element_text(size=10, face = "bold"),
      axis.title.y = element_text(size=12, angle=90, face = "bold", margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.text.x = element_text(size=10, face = "bold", margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.title.x = element_text(size=12, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, 'points'),
      strip.text.x = element_text(size = 5, color = "black", face = "bold.italic"),
      legend.margin = margin(0, 0, 0, 0), 
      legend.box.spacing = unit(0, "pt")
      # strip.background = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid")
    )
}


feature_list = colnames(rawdata[3:length(rawdata)])
rawdata$`Timepoint (day)`

## New column to group timepoints
values <- c("pre", "post", "post", "post")
index <- c(0, 12, 13, 14)
rawdata$Timepoint <- values[match(rawdata$`Timepoint (day)`, index)]


run_wilcoxon_summary <- function(dataframe, time_col, feature_cols, paired = TRUE) {
  
  ## Time point can be either pre/post for paired data or other grouping (treatment group?)
  ## paired set to true, option to include 'paired = FALSE' as a fourth argument for the function
  
  pvalues <- c()
  med_pre <- c()
  med_post <- c()
  
  for (feature in feature_cols) {
    # Ensure numeric
    dataframe[[feature]] <- as.numeric(dataframe[[feature]])
    
    # Get medians
    medians <- aggregate(dataframe[[feature]] ~ dataframe[[time_col]], dataframe, median)
    medians <- setNames(medians[,2], medians[,1])
    
    med_pre <- c(med_pre, medians["pre"])
    med_post <- c(med_post, medians["post"])
    
    # Wilcoxon test
    p <- wilcox.test(dataframe[[feature]] ~ dataframe$Timepoint, paired = paired)$p.value
    pvalues <- c(pvalues, p)
  }
  
  return(data.frame(
    feature = feature_cols,
    p_value = pvalues,
    median_pre = med_pre,
    median_post = med_post
  ))
}
results <- run_wilcoxon_summary(rawdata, "Timepoint", feature_list)

write.csv(results, "~/Path/To/Output/.csv", row.names=FALSE)







create_boxplots <- function(dataframe, feature_cols, group_col, paired = TRUE, custom_theme_fn = custom_theme) {
  plots_list <- list()
  
  my_comparisons <- list(c("pre", "post")) ## Can change to what is in the group_col
  
  for (feature in feature_cols) {
    dataframe[[feature]] <- as.numeric(dataframe[[feature]])
    
    plot <- ggplot(dataframe, aes(x = factor(.data[[group_col]]), y = .data[[feature]])) +
      geom_boxplot() +
      custom_theme_fn() +
      xlab("") +
      ylab("units") +
      ggtitle(feature) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=TRUE, size=3,
                                              step.increase = .075, tip.length = 0.01, hide.ns = TRUE) +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.6)
    
    plots_list[[feature]] <- plot
  }
  
  return(plots_list)
}

to_plot <- create_boxplots(rawdata, feature_list, "Timepoint")




save_plots_pdf <- function(plots_list, output_path) {
  pdf(output_path)
  
  plot_num <- 1
  plot_names <- names(plots_list)
  
  while (plot_num <= length(plots_list)) {
    plots_subset <- plots_list[plot_names[plot_num:min(plot_num + 3, length(plots_list))]]
    print(do.call(aplot::plot_list, c(unname(plots_subset), nrow = 2, ncol = 2)))
    plot_num <- plot_num + 4
    print(plot_num)
  }
  
  dev.off()
}

save_plots_pdf(to_plot, "~/Path/To/Output/BoxPlots.pdf")




















                      