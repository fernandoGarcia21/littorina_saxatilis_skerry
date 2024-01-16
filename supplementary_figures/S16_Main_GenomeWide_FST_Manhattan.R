################# Genome-wide FST Manhattan by Year ############################
#Generates a Manhattan plot for a specific year with the FST values of the 
#Skerry vs merged Crab ecotype population including both collinear loci and inversions. 
#Multiple years can be specified, so the result will be a grid of Manhattan plots.
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################

#Set the working directory
setwd('')

#Load the required libraries and packages
library(ggplot2)
library(dplyr)
require(scales)
library(hrbrthemes)
library("cowplot")
library(stringr)
library(qqman)
library(ggrepel)
library(forcats)
source('../generals/GeneralsSkerry.R')

#All datasets and functions are loaded from this package from Figure_S2:
source('../figure_2/GenomeWide_FST_Manhattan_by_Year.R') 

#Manhattan plot of the rest of the years for supplementary material
plt_1992 = generate_fst_manhattan_year('1992', "DO.92", "C 1992")
plt_2005 = generate_fst_manhattan_year('2005', "SKR.05", "S 2005")
plt_2018 = generate_fst_manhattan_year('2018', "ExpSkerry", "S 2018")

#Create a PDF file with the output plot
pdf('Manhattan_1992_2005_2018_Final_Supplementary_2024.pdf', width = 10, height=10)

#Create a grid with one plot on each row
prow <- plot_grid(
  plt_1992 + theme(legend.position="none"),
  plt_2005 + theme(legend.position="none"),
  plt_2018 + theme(legend.position="none"),
  nrow = 3,
  ncol = 1
)

# extract a legend that is laid out horizontally
legend_b <- get_legend(
  plt_1992 + 
    theme(legend.position = "bottom")
)
plot_grid(
  prow,
  legend_b,
  ncol = 1,
  rel_heights = c(1, .1))
dev.off()
