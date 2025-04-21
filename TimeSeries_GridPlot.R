library(ggplot2)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(tidyverse)
library(readxl)
library(viridis)   
library(RColorBrewer)
library(ggbreak)
library(cowplot)
library(patchwork)

## In case any variables are saved previously in environment you can run:
# rm(list=ls())
custom_theme <- function() {
  theme_classic() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size=9, face = "bold"),
      #legend.position = "none",
      legend.title = element_text(size=8),
      legend.text  = element_text(size = 6),
      legend.key.spacing.y = unit(0.0005, "cm"),
      legend.key.size = unit(0.4, "cm"),
      axis.text.y = element_text(size=6, face = "bold"),
      axis.title.y = element_text(size=8, angle=90, face = "bold"),
      axis.text.x = element_text(size=6, angle = 75, face = "bold", margin = margin(t = 0, r = 6, b = 0, l = 0)),
      axis.title.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(0, 'points'),
      strip.text.x = element_text(size = 5, color = "black", face = "bold.italic"),
      # strip.background = element_rect(color="black", fill="#FC4E07", size=1.5, linetype="solid")
      axis.title.x.top = element_blank(),
      axis.text.x.top  = element_blank(),
      axis.ticks.x.top = element_blank()
      
    )
}



## load data
Load_data <- read_excel("/Path/To/Input/Data/.xlsx", sheet = "%PBMC+0.0001")
## change 'Baseline' to 0
Load_data$time
Load_data$response


## Row 1 = headers
colnames(Load_data)
## Raw (y axis = % PBMC)
Raw <- Load_data[1:44,]


## Remove the NA for Raw$time?
Raw <- Raw[!is.na(Raw$time), ]
Raw$time
Raw$response


column_list = names(Load_data[6:length(Load_data)]) ## list of all columns that we want to make a plot for



###########################
##### Set parameters ######
###########################
biomarker_list = column_list
plot_list = list()

## Groups = PT
## X axis = 'time' (with split between 8 and 15 weeks)
## Y axis = 'biomarker'


for (item in biomarker_list) {

  Raw[[item]] = as.numeric(Raw[[item]])
  
  Days_of_interest = c(1, 15, 43, 56)
  Raw = Raw[(Raw$time %in% Days_of_interest),]
  
  tick_labels = sort(as.numeric(unique(Raw$time)))
  
  ### Time read in as.numeric not as qualitative
  plot <- ggplot(Raw, aes(x = as.numeric(time), y=as.numeric(.data[[item]]), color= factor(PT))) + 
    geom_point(aes(color=factor(PT), shape = factor(response)))  + ggtitle(paste(item)) +
    custom_theme() + ylab("% of PBMC") + xlab("weeks") + # scale_color_viridis(discrete=TRUE) +
    geom_line(aes(group = factor(PT)), linewidth = .4) + 
    # scale_x_break(breaks = c(18, 40), scales=.3, ticklabels= c(18, 40)) +
    scale_x_continuous(breaks=tick_labels, limits = c(0, max(tick_labels) + 2)) +
    scale_x_break(breaks = c(15, 42), scales=.3, ticklabels = tick_labels) + ##ticklabels only after the break
    labs(color = "PT Group", shape = "Response Status") + 
    scale_shape_manual(values = c(1, 2, 3, 4, 5)) + 
    guides(color = guide_legend(ncol=2))
  
  # plot
  
  pltName <- paste( 'Raw_plot_of_', item, sep = '' )
  plot_list[[pltName]] = plot
  

}


### Four per page
# initialize plotting device 
pdf("/Path/To/Output/PDF/.pdf")

plot_num <- 1

while (plot_num <= length(plot_list)) {
  
  remaining <- length(plot_list) - plot_num + 1
  
  if (remaining >= 4) {
    print(aplot::plot_list(plot_list[[plot_num]], plot_list[[plot_num + 1]], plot_list[[plot_num + 2]], plot_list[[plot_num + 3]], nrow = 2, ncol = 2))
    plot_num <- plot_num + 4
    
  } else if (remaining == 3) {
    print(aplot::plot_list( plot_list[[plot_num]], plot_list[[plot_num + 1]], plot_list[[plot_num + 2]], nrow = 2, ncol = 2))
    plot_num <- plot_num + 3
    
  } else if (remaining == 2) {
    # grid.arrange()
    print(aplot::plot_list(plot_list[[plot_num]], plot_list[[plot_num + 1]], nrow = 2, ncol = 2))
    plot_num <- plot_num + 2
    
  } else {
    print(aplot::plot_list(plot_list[[plot_num]], nrow = 2, ncol = 2))
    plot_num <- plot_num + 1
  }
}

dev.off()



