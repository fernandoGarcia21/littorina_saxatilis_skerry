############################## Figure 1  ###############################
# Figure 1: Trajectories of phenotypes and FST of collinear loci
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
source('./ShellLengthTrajectory.R')
source('./Phenotypes_Trajectories.R')
source('./ShellShapeScatterTrajectory.R')
source('./FSTTrayectories.R')

line_width = 1.0 #Width of the lines in the line plots of the trajectories
point_size = 1.5 #Size of the points in the trajectories
show_legend_geneal = TRUE #Either show or hide the legends on each trajectory plot

#Generate the lists of ggplot objects with the trajectories of the phenotypes and FST for Figure 1
list_shell_length_plt = generate_shell_length_plots() #Shell length scatter plot and error bars
list_phenotypes_plt = generate_phenotypes_plots() #Line plot with trajectory of four phenotypic traits
list_shell_scatter_plt = generate_phenotypes_scatter_plots() #Scatter plot based shell length and ln(gw) of the shell
tmp_plots_scatter_skerry = list_shell_scatter_plt[c('1992','1996','2002','2005','2018','2021')] #set the order of the years for the trajectory
list_FST_plt = generate_FST_plots() #FST trajectory of collinear loci

#Generate the list of ggplot objects with the trajectories of the inversion arrangements for Figure 3
list_inversions_plt = generate_inversions_plots()


##################### Figure 1 #################################################
#Generate a PDF with the Figure 1 trajectories
pdf('Original Fig1_Trajectories_combined_Final_elongated.pdf', width = 5, height = 10)
#Plot a grid with two columns: 1-Skerry trajectories. 2-Wave trajectories
plot_grid(list_shell_length_plt[['Skerry']], list_shell_length_plt[['Wave']],
          list_phenotypes_plt[[1]][['Skerry']], list_phenotypes_plt[[1]][['Wave']],
          list_phenotypes_plt[[2]][['Skerry']], list_phenotypes_plt[[2]][['Wave']],
          plot_grid(plotlist = tmp_plots_scatter_skerry, nrow = 1), list_shell_scatter_plt[['Wave']],
          list_FST_plt[['Skerry']], list_FST_plt[['Ref']],
          ncol = 2, nrow = 5, rel_widths = c(0.85,0.15),
          rel_heights = c(0.24,0.23,0.23,0.05,0.25))
dev.off()
################################################################################
