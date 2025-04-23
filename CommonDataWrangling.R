## load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

## Tidyverse and base R solutions
dataframe <- read.csv('/Users/lynchmt/Documents/OXEL/OXELdatabase_Update.01.03.24.csv')
column <- colnames(dataframe)[75] 
print(colnames(dataframe)[1:100])
## Use dataframe[[column]] when column is stored as a string variable: column <- "column_name"
## Otherwise, dataframe$column_name

########################################################################
############################### Cleaning ###############################
########################################################################

## Missing Data
dataframe_temp <- dataframe %>% drop_na() #tidyverse
dataframe_temp <- dataframe[complete.cases(dataframe), ] # base R

## Remove row when a specific column has missing data
dataframe_temp <- dataframe %>% filter(!is.na(.data[[column]])) # tidyverse
dataframe_temp <- dataframe[!is.na(dataframe[[column]]), ] # base R


## Remove white space
dataframe_temp <- dataframe %>% mutate(across(everything(), str_trim)) #tidyverse
dataframe_temp[[column]] <- trimws(dataframe[[column]]) # base R


## Format column names
# All lower case
dataframe_temp <- dataframe %>% rename_with(tolower) #tidyverse
names(dataframe_temp) <- tolower(names(dataframe)) # base R

# All upper case
names(dataframe_temp) <- toupper(names(dataframe))

# Remove character (e.g., remove ".x")
names(dataframe_temp) <- gsub("\\.x", "", names(dataframe))

# Add suffix
names(dataframe_temp) <- paste0(names(dataframe), "_suffix")


##############################################################################
############################### Pre-Processing ###############################
##############################################################################

## Filter
dataframe_temp <- dataframe %>%  filter(Arm.assigned == "A Nivo") #tidyverse
dataframe_temp <- dataframe[dataframe$Arm.assigned == "A Nivo", ] # base R

## Sorting
dataframe_temp <- dataframe %>% arrange(Arm.assigned) #tidyverse
dataframe_temp <- dataframe[order(dataframe[[column]]), ] # base R


example_df <- data.frame(
  ID = 1:4,
  Priority = c("Medium", "Low", "High", "Medium"),
  Score_Week1 = c(10, 15, 12, 13),
  Score_Week2 = c(14, 18, 13, 17),
  Score_Week3 = c(16, 20, 14, 19)
)

print(example_df)


## Sorting with factor levels


# Tidyverse:
sort_column <- "Priority"

example_temp <- example_df %>%
  mutate({{ sort_column }} := factor(.data[[sort_column]], levels = c("High", "Medium", "Low"))) %>%
  arrange(.data[[sort_column]])

print(example_temp)


# base R
sort_column <- "Priority"
# Convert to factor with custom order
example_df[[sort_column]] <- factor(example_df[[sort_column]], levels = c("Low", "Medium", "High"))
# Sort
example_temp <- example_df[order(example_df[[sort_column]]), ]
print(example_temp)


## Combine two datasets
# df_combined <- rbind(df1, df2) # base R
# df_combined <- bind_rows(df1, df2) # tidyverse

## Combine and summarize
aggregate(CCL4_12w ~ Arm.assigned, data = dataframe, FUN = mean) # base R






###############################################################
########################## Reshaping ##########################
###############################################################


## Wide to long
# Tidyverse
example_long <- example_df %>%
  pivot_longer(
    cols = starts_with("Score"),
    names_to = "Week",
    values_to = "Score"
  )

print(example_long)


# Base R
example_long <- reshape(
  example_df,
  varying = list(3:5),
  v.names = "Score",
  timevar = "Week",
  times = names(example_df)[3:5],
  direction = "long"
)

print(example_long)


## Long to wide
# Tidyverse
example_wide <- example_long %>%
  pivot_wider(names_from = Week, values_from = Score)

print(example_wide)


# Base R
example_wide <- reshape(
  example_long,
  idvar = "ID", timevar = "Week",
  direction = "wide"
)

print(example_wide)



##################################################################
########################### Formatting ###########################
##################################################################

## Rename columns
## OLD name = OS and NEW name = OVerallSurvival
dataframe_temp <- dataframe %>% rename(OverallSurvival = OS) #tidyverse

## OLD name = OS and NEW name = OVerallSurvival
names(dataframe)[names(dataframe) == "OS"] <- "OverallSurvival" # base R


## Re-order columns

# df <- df[c("col3", "col1", "col2")] # base R









