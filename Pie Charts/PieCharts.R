rm(list = ls())
library(geomtextpath)
library(ggplot2)
library(stringr)
library(dplyr)
library(webr)
library(RColorBrewer)
library(ggforce)
library(cowplot)
library(ggtext)
library(grid)
dataframe <- read.table(file = './PieChart/filename.tsv', sep = '\t', header = TRUE)

UMIcolumns <- grep("uniqueUMICountAggregated", colnames(dataframe), value = TRUE)
pt5 <- grep("pt5", UMIcolumns, value = TRUE)
pt6 <- grep("pt6", UMIcolumns, value = TRUE)
pt7 <- grep("pt7", UMIcolumns, value = TRUE)
pt8 <- grep("pt8", UMIcolumns, value = TRUE)
pt11 <- grep("pt11", UMIcolumns, value = TRUE)



## 

pie_theme <- function() {
  theme_void() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=10, face="bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
      legend.text = element_text(size = 6), # Adjusts the legend text size
      legend.position = "bottom",
      legend.key.size = unit(0.2, 'cm') ,    # Adjusts the size of the legend keys (the boxes next to the legend text)
      legend.margin = margin(t = 0, r = 0, b = 10, l = 0), # Adjusts the margin around the legend
      legend.title = element_blank(),
      panel.background = element_blank(), plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

pie_theme_2 <- function() {
  theme_void() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=10, face="bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
      legend.text = element_text(size = 8), # Adjusts the legend text size
      legend.position = "bottom",
      legend.key.size = unit(0.3, 'cm') ,    # Adjusts the size of the legend keys (the boxes next to the legend text)
      legend.margin = margin(t = 0, r = 30, b = 10, l = 30), # Adjusts the margin around the legend
      legend.title = element_blank(),
      panel.background = element_blank(), plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}


####################################################################################################
########################### Charts w/ top 5 highest frequency by patient ###########################
####################################################################################################


top5_function <- function(patient_list){
  return_list <- c()
  for (p_t in patient_list){
    data_temp <- data.frame(vGene_aa = paste(sapply(str_split(dataframe$vGene, "\\*"), `[`, 1), dataframe$aaSeqCDR3, sep="_"))
    data_temp[[p_t]] <- dataframe[[p_t]]
    data_temp$percentage <- dataframe[[p_t]] / sum(dataframe[[p_t]] )
  
    data_temp <- data_temp[order(data_temp$percentage , decreasing = TRUE),]  
  
    running_total = 0
    i = 1
    pie_chart_df = data_temp[0,] ## add the column names to pie_chart_df
  
    while(running_total < .25){
    pie_chart_df <- data_temp[0:i,]
    running_total = sum(pie_chart_df$percentage)
    # print(running_total)
    i = i + 1
  }
  
    ## Combine by vGene_aa
    value_cols <- c(p_t, 'percentage')
  
    top5 <- pie_chart_df$vGene_aa[1:5]
    print(top5)
    return_list <- c(return_list, top5)
  }
  return(return_list)
}
top20 <- top5_function(pt11)

### Color unique top 20
print(unique(top20))

start_palette <- scales::hue_pal()(length(unique(top20)))
names(start_palette) <- unique(top20)







#############################################################################################
########################### Pie chart colored by V gene, no Donut ###########################
#############################################################################################
create_pie_chart <- function(patient_list, start_palette){
  return_list <- c()
  for (p_t in patient_list){
    data_temp <- data.frame(vGene_aa = paste(sapply(str_split(dataframe$vGene, "\\*"), `[`, 1), dataframe$aaSeqCDR3, sep="_"))
    data_temp[[p_t]] <- dataframe[[p_t]]
    data_temp$percentage <- dataframe[[p_t]] / sum(dataframe[[p_t]] )
    
    data_temp <- data_temp[order(data_temp$percentage , decreasing = TRUE),]  
    
    running_total = 0
    i = 1
    pie_chart_df = data_temp[0,] ## add the column names to pie_chart_df
    
    while(running_total < .25){
      pie_chart_df <- data_temp[0:i,]
      running_total = sum(pie_chart_df$percentage)
      # print(running_total)
      i = i + 1
    }
    
    ## Combine by vGene_aa
    value_cols <- c(p_t, 'percentage')
    top5 <- pie_chart_df$vGene_aa[1:5]

    pie_chart_df <- pie_chart_df[order(pie_chart_df$percentage, decreasing = TRUE),]  
    annotate = pie_chart_df$percentage
    final_row <- data.frame(vGene_aa = "Other")
    final_row[[p_t]] <- 0
    final_row$percentage <- 1 - sum(pie_chart_df$percentage)
    pie_chart_df <- rbind(pie_chart_df, final_row)
    
    ## re-order:
    pie_chart_df$vGene_aa <- factor(pie_chart_df$vGene_aa, levels = rev(pie_chart_df$vGene_aa))
    
    
    # Create a column for legend variable
    pie_chart_df <- pie_chart_df %>%
      mutate(legend_category = ifelse(vGene_aa %in% top5, as.character(vGene_aa), NA))
    
    print(head(pie_chart_df))
    
    new_palette <- ifelse(names(start_palette) %in% pie_chart_df$legend_category,
                          start_palette,
                          "grey50")
    names(new_palette) <- names(start_palette)
    
    
    missing_colors <- setdiff(pie_chart_df$vGene_aa, names(new_palette))
    grey_colors <- rep("grey50", length(missing_colors))
    names(grey_colors) <- missing_colors
    
    final_palette <- c(new_palette, grey_colors)
    
    
    p <- ggplot(pie_chart_df, aes(x = "", y = percentage, fill = vGene_aa)) +
      geom_bar(stat = "identity", width = .001, color="white") +
      coord_polar(theta = "y") +
      labs(title = paste("Pie Chart with Top 5 in Legend\n", p_t)) +
      scale_fill_manual(values = final_palette,
                        breaks = top5) +
      pie_theme_2() +
      theme(legend.title = element_blank()) + guides(fill = guide_legend(nrow = 3)) 
    
    p <- p + 
      annotate("text",
        x = 1,        #  0.75 to get near top-right
        y = 1,      # radius
        label = paste(length(annotate)),
        size = 6.5,        # text size
        hjust = -1,       # left align horizontally
        vjust = 0, fontface = "bold")
    
    png(paste("./PieChart/ByPatient/", p_t,".png", sep=""), 
        width=4.25, height=5.5, units="in", res=300)
    
    print(p)
    dev.off()
    
    
  }
  return(return_list)
}
create_pie_chart(pt11, start_palette)

                                                                                                                                            











#################################################################################################

########################################## Pie-Donut Chart ######################################

#################################################################################################

pie_theme_3 <- function() {
  theme_void() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=5, face="bold", margin = margin(t = 15, r = 0, b = 0, l = 0)),
      legend.text = element_text(size = 3, margin = margin(b=2)), 
      legend.position = "right",
      legend.key.size = unit(0.05, 'cm') ,    # Adjusts the size of the legend keys (the boxes next to the legend text)
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # Adjusts the margin around the legend
      legend.spacing = unit(0.5, 'cm'),
      legend.title = element_blank(),
      panel.background = element_blank(), plot.background = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

#################################################################################################

####################################### Defining the palettes ###################################

#################################################################################################


### Assign colors outside of the loop so they stay consistent 
data_temp <- data.frame(vGene_aa = sapply(str_split(dataframe$vGene, "\\*|-"), `[`, 1))
data_temp[[UMIcolumns[1]]] <- dataframe[[UMIcolumns[1]]]
pie_chart_df <- aggregate(data_temp[, UMIcolumns[1]], by = list(data_temp[['vGene_aa']]), FUN = sum)
pie_chart_df <- pie_chart_df[order(pie_chart_df$x, decreasing = TRUE),]  

color_palette <- scales::hue_pal()(length(pie_chart_df$Group.1))
names(color_palette) <- pie_chart_df$Group.1

suffix_numbers <- as.numeric(sub(".*TRBV", "", names(color_palette)))
ordered_names_all <- names(color_palette)[order(suffix_numbers)]
color_palette <- color_palette[ordered_names_all]


### Choose color for TRBV6, TRBV10:
color_palette[['TRBV6']] = "#FB61D7"
color_palette[['TRBV17']] = "#C49A00"
color_palette[['TRBV10']] = "#7CAE00"

brown_white_palette <- colorRampPalette(c(color_palette[['TRBV6']], "white"))
green_white_palette <- colorRampPalette(c(color_palette[['TRBV10']], "white"))

### OUTER == NA for all except V6, V10 ??
data_outer <- data.frame(inner = sapply(str_split(dataframe$vGene, "\\*|-"), `[`, 1),
                         outer = sapply(str_split(dataframe$vGene, "\\*"), `[`, 1))

data_outer <- data_outer[data_outer$inner %in% c('TRBV6', 'TRBV10'), ]
data_outer <- unique(data_outer)

# Split by inner_name
outer_split <- split(data_outer, data_outer$inner)
V10_colors <- green_white_palette(length(outer_split$TRBV10$outer)+1)[1:length(outer_split$TRBV10$outer)]
names(V10_colors) <- outer_split$TRBV10$outer
V6_colors <- brown_white_palette(length(outer_split$TRBV6$outer)+1)[1:length(outer_split$TRBV6$outer)]
names(V6_colors) <- outer_split$TRBV6$outer


suffix_numbers <- as.numeric(sub(".*TRBV6-", "", names(V6_colors)))
ordered_names_V6 <- names(V6_colors)[order(suffix_numbers)]
V6_colors <- V6_colors[ordered_names_V6]


suffix_numbers <- as.numeric(sub(".*TRBV10-", "", names(V10_colors)))
ordered_names_V10 <- names(V10_colors)[order(suffix_numbers)]
V10_colors <- V10_colors[ordered_names_V10]



#################################################################################################
########################################## Loop by patient #####################################
#################################################################################################


patient = pt5
png("./PieChart/Plot_Patient5.png", width = 400 * (length(patient) + 1), height = 400, res=400)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(patient) + 1)))

i = 1
for (p_t in patient){
  
  data_inner <- data.frame(inner = sapply(str_split(dataframe$vGene, "\\*|-"), `[`, 1))
  head(data_inner)
  data_inner[[p_t]] <- dataframe[[p_t]]
  print(head(data_inner))
  
  data_inner <- aggregate(data_inner[, p_t], by = list(data_inner[['inner']]), FUN = sum)
  names(data_inner) <- c('inner_name', 'inner_count')
  data_inner <- data_inner[order(data_inner[['inner_count']], decreasing = TRUE),]  
  data_inner <- data_inner[(data_inner[['inner_count']] != 0),]  
  
  ## Calculate fraction
  data_inner$fraction <- data_inner$inner_count / sum(data_inner$inner_count)
  ## Calculate ymax
  data_inner$ymax <- cumsum(data_inner$fraction)
  data_inner$ymin <- c(0, head(data_inner$ymax, -1))
  
  data_inner$inner_name <- factor(data_inner$inner_name, levels = rev(data_inner$inner_name))
  
  label_radius <- 1.05
  
  ### coordinates for labels
  data_inner <- data_inner %>%
    mutate(
      mid_fraction = (ymin + ymax) / 2,
      theta = pi/2 - 2 * pi * mid_fraction,
      label_x = .8 * cos(theta),
      label_y = .8 * sin(theta),
     #  hjust = ifelse(cos(theta) > 0, .7, .3),
      # vjust = ifelse(sin(theta) > 0, 0, 1),
      perc = scales::percent(fraction / sum(fraction), accuracy = 0.1)
    )
  
  

  ### OUTER == NA for all except V6, V10 ??
  data_outer <- data.frame(inner = sapply(str_split(dataframe$vGene, "\\*|-"), `[`, 1),
                           outer = sapply(str_split(dataframe$vGene, "\\*"), `[`, 1))
  data_outer[[p_t]] <- dataframe[[p_t]]
  head(data_outer)
  
  data_outer <- data_outer[data_outer$inner %in% c('TRBV6', 'TRBV10'), ]
  
  data_outer <- aggregate(data_outer[, p_t], by = list(data_outer[['outer']], data_outer[['inner']]), FUN = sum)
  head(data_outer)
  names(data_outer) <- c('outer_name', 'inner_name', 'outer_count')
  
  
  data_outer <- merge(data_outer, data_inner[, c("inner_name", "ymin", "ymax")], 
                      by = "inner_name", all.x = TRUE)
  
  # Split by inner_name
  outer_split <- split(data_outer, data_outer$inner_name)
  
  # Loop over each group
  outer_split <- lapply(outer_split, function(df) {
    df$sub_fraction <- df$outer_count / sum(df$outer_count)
    df$sub_ymax <- df$ymin[1] + 
      (df$ymax[1] - df$ymin[1]) * cumsum(df$sub_fraction)
    df$sub_ymin <- df$ymin[1] + 
      (df$ymax[1] - df$ymin[1]) * c(0, head(cumsum(df$sub_fraction), -1))
    df
  })
  
  # Recombine into a single dataframe
  data_outer <- do.call(rbind, outer_split)
  
  data_labels <- rbind(data_inner[1:5,], data_inner[data_inner$inner_name %in% c('TRBV6', 'TRBV10'), ])
  
  
  p <- ggplot() +
    # Inner pie
    geom_arc_bar(data = data_inner,
      aes(x0 = 0, y0 = 0,
        r0 = 0, r = 1,
        start = 2 * pi * ymin,
        end = 2 * pi * ymax,
        fill = inner_name),
      color = "white") + coord_fixed()  + # Outer ring
    geom_arc_bar(data = data_outer,
      aes(x0 = 0, y0 = 0,
        r0 = 1, r = 1.2,
        start = 2 * pi * sub_ymin,
        end = 2 * pi * sub_ymax,
        fill = outer_name),
      color = "white") +
    scale_fill_manual(values = c(color_palette, V6_colors, V10_colors), breaks = c(ordered_names_all, ordered_names_V6, ordered_names_V10),
                      drop = FALSE) +
    coord_fixed() +
    pie_theme_3() +
    labs(title = paste(str_split(p_t, "_")[[1]][1])) +
    guides(fill = guide_legend(
      order = 1,
      ncol = 3
    ))  +
    geom_textpath(data = data_labels,
                 # hjust = data_labels$hjust,
                 # vjust = data_labels$vjust,
                  aes(x = label_x, y = label_y, 
                      label = paste0(inner_name, "\n (", perc, ")")),
                  size = .75, show.legend = FALSE, fontface = "bold"
    )
  
  
  legend <- get_legend(p + theme(legend.position = "right"))
  p <- p + theme(legend.position = "none")
  print(p, vp = viewport(layout.pos.row = 1, layout.pos.col = i))
  i = i + 1
}

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = length(patient) + 1))
grid.draw(legend)
popViewport()
dev.off()













