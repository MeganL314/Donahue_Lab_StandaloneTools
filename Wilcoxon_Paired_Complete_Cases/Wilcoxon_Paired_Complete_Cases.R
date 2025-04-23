library(tidyverse)
library(dplyr)
# import raw data
# to create this, use the import feature in R studio

library(readxl)
rawdata <- read_excel("~/Path/To/Excel/File.xlsx", 
                      sheet = "All data", col_types = c("text", 
                                                        "text", "text", "text", "text", "text", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric", 
                                                        "numeric", "numeric", "numeric"))
View(rawdata)
# pairs of variable comparisons. THE ORDER OF THE PAIRS MUST MATCH THE ORDER OF THE NAMES ABOVE
pairs <- list(
  c(1,15),
  c(1,29),
  c(1,43)
)

# start and end column index
start_col <- 7
end_col <- 15


Patient_ID_Column = "Patient #"


#!!! probably do not need to change anything below!!!


# create an empty results data frame to append each experiment to
# results <- data.frame()
# for loop that loops over all experiments

wilcox_pvalue_list = c()
analyte_list = c()
group_list = c()
sample_size_list = c()
median1_list = c()
median2_list = c()
results = data.frame()
for(pair in pairs){
  # extract the name and associated number of the variable
  num1 <- pair[1]
  num2 <- pair[2]
  # subset the raw data to the two groups for the current experiment
  temp <- subset(rawdata, Group == num1 | Group == num2)
  temp <- as.data.frame(temp)
  
  median1 <- sapply(subset(temp, Group == num1)[ ,start_col:end_col], FUN = median, na.rm = TRUE)
  median1_list <- c(median1_list, median1)
  
  median2 <- sapply(subset(temp, Group == num2)[ ,start_col:end_col], FUN = median, na.rm = TRUE)
  median2_list <- c(median2_list, median2)
  
  
  for (col in start_col:end_col){
    print(col)
    list_rownames <- c()
    i = 1
    for (i in 1:length(subset(temp, Group == num1)[[1]])){
      first_value = temp[,col][i]
      first_id = temp[,1][i]
      # print(first_value)
      pair = i + length(subset(temp, Group == num1)[[1]])
      second_value = temp[,col][pair]
      second_id = temp[,1][pair]
      if (is.na(first_value) | is.na(second_value)){
        # remove rows from temp
        list_rownames = c(list_rownames, first_id[[1]], second_id[[1]])
        # temp <- subset(temp, 'Sample:' %in% c(first_id[[1]], second_id[[1]]))
      }}
    print(list_rownames)
    
    
    new_temp <- temp[ !temp[[Patient_ID_Column]] %in% list_rownames, ]
    print(length(row.names(new_temp)))
    
    if (length(new_temp[,col]) >= 2){
      wilcox_temp <-  wilcox.test(new_temp[[col]] ~ as.factor(new_temp$Group), alt="two.sided", paired=T, correct=T, conf.int=F, conf.level=0.95)$p.value
      
      wilcox_pvalue_list = c(wilcox_pvalue_list, wilcox_temp)
      analyte_list = c(analyte_list, colnames(new_temp[col]))
      group_list = c(group_list, paste("Groups ", num1, " and ", num2, sep = ""))
      sample_size_list = c(sample_size_list, paste("n= ", length(new_temp[,col]), sep=""))
    }
    else{
      wilcox_pvalue_list = c(wilcox_pvalue_list, "No value")
      analyte_list = c(analyte_list, colnames(new_temp[col]))
      group_list = c(group_list, paste("Groups ", num1, " and ", num2, sep = ""))
      sample_size_list = c(sample_size_list, paste("n= 0"))
      
    }}
  # i = i + 1
  #wilcox_results <- data.frame(wilcox_pvalue_list, analyte_list, group_list, sample_size_list)
  #tmp_results <- cbind(wilcox_results, median1, median2)
  #results <- rbind(tmp_results, results)
}

results <- data.frame(wilcox_pvalue_list, analyte_list, group_list, sample_size_list, median1_list, median2_list)

print(results)

write.csv(results, "~/path/to/output.csv", row.names=FALSE)
