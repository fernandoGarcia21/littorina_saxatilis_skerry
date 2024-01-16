######################## Figure 2: Pies charts and Manhattan Plot#########
# Figure 2A: Pie charts of collinear loci outside the neutral expectation
# Figure 2D: Manhattan plot of FST for collinear loci and inversions 2021
#
#author: "Diego Garcia"
#date: "2023-05-11"
#########################################################################

#Set the working directory
setwd('.')

#Load libraries
library(gridExtra)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(cowplot)
library(hrbrthemes)
library(stringr)
library(qqman)

#Load R scripts with the necessary functions to generate the different
# analysis and the corresponding figures
source('../generals/GeneralsSkerry.R')
source('PieChartsEnvelope.R')
source('GenomeWide_FST_Manhattan_by_Year.R')

#Invoke the function that generates the pies chart of the neutral expectation
grid_envelope_pie_charts = generate_envelope_pie_charts(FALSE)
#Invoke the function that generates the Manhattan plot of FST for skerry 2021
fst_manhattan_2021 = generate_fst_manhattan_year('Genome-wide FST in Skerry 2021 vs. the average Crab ecotype', "SRK.21", "S 2021")

#Create a PDF file with the output figure
pdf('Figure2_inversions_imputed_newpercentages.pdf', width = 10, height = 5)

#Create a grid with two rows: Pie charts and Manhattan plot respectively
plot_grid(grid_envelope_pie_charts, 
          fst_manhattan_2021,
          ncol=1,
          nrow=2,
          rel_heights = c(0.4,0.6))

dev.off()


