######################## Figure 1 & Figure 3 ############################
# Figure 3: Trajectories of the inversion arrangements
#
#author: "Diego Garcia"
#date: "2023-05-11"
#########################################################################

#Set the working directory
setwd('.')

#Load libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)
library(rlang)

#Load R scripts with the necessary functions to generate the different
# analysis and the corresponding figures
source('../generals/GeneralsSkerry.R')
source('./AllInversionsTrajectoriesAndDrift.R')

line_width = 1.0 #Width of the lines in the line plots of the trajectories
point_size = 1.5 #Size of the points in the trajectories
show_legend_geneal = TRUE #Either show or hide the legends on each trajectory plot

#Generate the list of ggplot objects with the trajectories of the inversion arrangements for Figure 3
list_inversions_plt = generate_inversions_plots()

##################### Figure 3 #################################################
#Generate a PDF with the Figure 3 inversion arrangement trajectories
pdf('LS_Fig3_InversionsTrajectories_elongated_imputed.pdf', width = 5, height = 10)
#Plot a grid with two columns: 1-Skerry trajectories. 2-Wave trajectories
plot_grid(list_inversions_plt[[1]][['Skerry']], list_inversions_plt[[1]][['Wave']],
          list_inversions_plt[[2]][['Skerry']], list_inversions_plt[[2]][['Wave']],
          list_inversions_plt[[3]][['Skerry']], list_inversions_plt[[3]][['Wave']],
          NULL,
          ncol = 2, nrow = 3, rel_widths = c(0.85,0.15), rel_heights = c(0.3,0.25,0.25))
dev.off()
################################################################################