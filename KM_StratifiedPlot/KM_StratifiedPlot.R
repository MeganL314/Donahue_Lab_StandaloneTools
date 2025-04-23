#install.packages("readxl")
#install.packages("dplyr")
#install.packages("ggplot")
#install.packages("ggplot2")
#install.packages("tidyverse")
#install.packages("tidyr")
#install.packages("survminer")
#install.packages("gtsummary")
#install.packages("flextable")

library(ggplot2)
library(tidyr)
library(rstatix)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(pROC)
library(nnet)
library(ROCR)
library(survival)
library(survminer)
library(tidyverse)

library(readxl)
library(tidyverse)
library(survival)
library(survminer)
# library(gtsummary)
# library(flextable)


### Change font size, header position, etc. 
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, face="bold", size=18),
      legend.key.height = unit(.4, 'cm'),
      legend.key.width = unit(.4, 'cm'),
      legend.text=element_text(size=15),
      legend.title=element_text(size=15),
      axis.text.y = element_text(size=16),
      axis.text.x = element_text(size=16),
      axis.title.y = element_text(size=17, angle = 90,  face="bold", margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.title.x = element_text(size=17, face="bold"),
      
    )
}

## load data
dataframe <- read_excel("data.xlsx", sheet = "sheet1")

## create list of analyte names using column names from dataframe
list = names(dataframe[10:length(dataframe)])


####################################################################################
#################################### Functions: ####################################
####################################################################################

create_low_high_list <- function(feature_data, cutoff) { 
  
  ## feature_data: a numeric vector of feature values (e.g., dataframe[[feature]])
  ## cutoff: a value between 0 and 1 indicating the percentile cutoff (e.g., 0.5 = median)

  threshold <- quantile(feature_data, probs = cutoff, na.rm = TRUE)
  low_high_list <- ifelse(feature_data >= threshold, 'Upper', 'Lower')
  return(list(labels = low_high_list, threshold = threshold))
  
}


create_labels <- function(data_for_labels, df_followup) { 
  label_list <- c()
  names_split <- strsplit(names(data_for_labels), "=")
  group_names <- sapply(names_split, function(x) x[2])
  data_for_labels =  unlist(data_for_labels)
  
  for (i in 1:length(data_for_labels)){
    if (is.na(data_for_labels[i])){
      label = paste(group_names[i],": Median = not reached at ", max(df_followup), " days", sep="")
    }else{label = paste(group_names[i],": Median = ", data_for_labels[i], " days", sep="")}
    label_list <- c(label_list, label)
  }
  return(label_list)
}






#################################################################################
################################ Set parameters #################################
#################################################################################
my_analytes = headers_character_unlist
dataframe = dataframe
outcome_time = "followup"
outcome_event = "dead"



create_KM_plots <- function(my_analytes, dataframe, outcome_time, outcome_event) { 
  
  p_val_list = c() ## initialize character vector to store p values
  median_high = c() ## initialize character vector to store median value for 'Upper' group
  median_low = c() ## initialize character vector to store median value for 'Lower' group
  feature_name = c() ## initialize character vector to store protein/analyte/biomarker/feature name
  
  
  for (analyte in my_analytes) {
    
    # dataframe_temp = dataframe %>% drop_na(analyte)
    dataframe_temp <- dataframe[!is.na(dataframe[[analyte]]), ]
    
    
    # call function to create labels of 'Upper' and 'Lower', also returns the cut-off label to use for KM plot
    result <- create_low_high_list(dataframe_temp[[analyte]], 0.5) ## This is where you specify 0.5 = Median, .75 = upper quantile, etc.
    low_high_list <- result$labels
    label_c <- paste("Cutoff:", round(as.numeric(result$threshold), digits=3), sep = " ")
  
    stratify_df = data.frame("time" = dataframe_temp[[outcome_time]], "event" = dataframe_temp[[outcome_event]], "group" = low_high_list)
    stratify_df$group <- factor(stratify_df$group, levels = c("Lower", "Upper"))
    
    fit <- surv_fit(Surv(time, event)~group, data = stratify_df)
    p_val <- surv_pvalue(fit, data=NULL, method = "survdiff", test.for.trend = FALSE, combine = TRUE)$pval
  
    data_for_labels  <- summary(fit)$table[c("group=Lower", "group=Upper"), "median"]
    legend_labels <- create_labels(data_for_labels, dataframe_temp[[outcome_time]])
  
    ### Update character vectors ###
    feature_name <- c(feature_name, analyte)
    p_val_list <- c(p_val_list, p_val)
    median_high <- c(median_high, summary(fit)$table["group=Upper", "median"])
    median_low <- c(median_low, summary(fit)$table["group=Lower", "median"])
  
    ### Plot stratified KM - 1 per page ###
    # Define the plot
    ggsurv <- ggsurvplot(fit, data = stratify_df, risk.table = TRUE, pval = TRUE, 
                         title=paste(as.character(analyte), " Plot",sep = ""),
                         ggtheme = custom_theme(), pval.size = 5, legend.labs = legend_labels, 
                         tables.y.text = FALSE, pval.coord = c(90, 0.2), palette = c("#ff9a90", "#00A170"))
    
    ggsurv$plot <- ggsurv$plot + guides(fill = guide_legend(nrow = 2), color = guide_legend(nrow = 2))
    
    ggsurv$plot <- ggsurv$plot + 
      ggplot2::annotate("text", x=90, y=.14, label = label_c, size = 5, hjust=0, fontface=1) ## you will need to adjust coordinates (x=90 and y=0.14)
  
  
    print(ggsurv)
  
  }
  return(list(
    feature_name = feature_name,
    p_val_list = p_val_list,
    median_high = median_high,
    median_low = median_low
  ))
  
  
  }



outfile_name = "ArmA"
pdf(paste("/Path/To/Save/Plot/KM_Test_", outfile_name, ".pdf", collapse = ""))

results_list <- create_KM_plots(my_analytes, dataframe, outcome_time, outcome_event)

dev.off()

results <- cbind(results_list$feature_name, results_list$p_val_list, results_list$median_high, results_list$median_low)


write.table(results, file = "/Path/To/Save/csv/.csv", sep = ",", row.names = F, quote=F)




