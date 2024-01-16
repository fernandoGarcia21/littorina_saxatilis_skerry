#################### Space-Time FST Correlation #################
# Figure 2B: FST correlation for Spatial Outliers
# Figure 2C: FST correlation for inversions
#
#author: "Diego Garcia"
#date: "2024-01-04"
#########################################################################

#Set the working directory
setwd('.')

#Load libraries
library(ggplot2)
library(dplyr)
require(scales)
library(hrbrthemes)
library("cowplot")
library('stringr')
library(plotly)

#Define the path and name of the files that contain FST values for space and time
#These files are created in the R script 2_seqsnptm009_time vs space.R
fst_correlations_file_outliers <- "../Data/SEQSNPTM009_TIME VS SPACE_OUT.txt" #collinear outliers
fst_correlations_file_inversions <- "../Data/SEQSNPTM009_TIME VS SPACE_INV.txt" #inversions

#Read file and load data
fst_correlations_outliers_data <- read.table(fst_correlations_file_outliers, header = T, check.names = FALSE)
fst_correlations_inversions_data <- read.table(fst_correlations_file_inversions, header = T, check.names = FALSE)

#Define the colour and shape scheme for the correlation plot
colors_fst = c('MoralesEtAl' = '#73BADA', 'WestramEtAl' = '#E09157', 'R' = '#7A659E', 'Inversion' = "#F34573")
shapes_fst = c('MoralesEtAl' = 16, 'WestramEtAl' = 17, 'Inversion' = 15)


##################### OUTLIERS CORRELATION PLOT ##############################################
#Linear model of correlation
m_out = lm(fst_skerry~fst_area, fst_correlations_outliers_data) #slope
r2_out = format(summary(m_out)$r.squared, digits = 2) #R2 or the coefficient of determination

#Estimate the correlation coefficient and p value
corr <- cor.test(x=fst_correlations_outliers_data$fst_area, y=fst_correlations_outliers_data$fst_skerry, method = 'spearman')
corr

#Generate a ggplot object with the space-time correlation for OUTLIERS (collinear loci)
plt_outliers = ggplot(fst_correlations_outliers_data, 
              aes(x=fst_area, y=fst_skerry))+
  geom_point(alpha=0.7, size=2, aes(color=factor(cat),
                                      shape=factor(cat))) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color=colors_fst['R'], linetype='F1', size=0.8)+
  annotate(geom="text",x=0.9,y=0.6,
           label=(paste0("R2==",r2_out)),
           parse=TRUE, color=colors_fst['R'], size=3, fontface = 'bold')+
  ylim(-0.1, 1)+
  xlim(0, 1)+
  xlab('FST in space (Crab-Wave in nearby location)')+
  ylab('FST in time (Skerry 1992 vs 2021)')+
  theme_minimal()+
  #ggtitle(p_title)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10),
        #plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
  scale_color_manual(values = colors_fst)+
  scale_shape_manual(values = shapes_fst)

#Generate a PDF file with the correlation plot
pdf(width = 5, height = 3, file = './FST_Space_Time/Space_Time_FST_Outliers_2024.pdf')
plt_outliers
dev.off()
########################################################################################################


############################# INVERSIONS CORRELATION PLOT ##############################################
#Linear model of correlation
m_inv = lm(fst_skerry~fst_area, fst_correlations_inversions_data) #slope
r2_inv = format(summary(m_inv)$r.squared, digits = 2) #R2 or the coefficient of determination

#Estimate the correlation coefficient and p value
corr <- cor.test(x=fst_correlations_inversions_data$fst_area, y=fst_correlations_inversions_data$fst_skerry, method = 'spearman')
corr

#Generate a ggplot object with the space-time correlation for OUTLIERS (collinear loci)
plt_inversions = ggplot(fst_correlations_inversions_data, 
             aes(x=fst_area, y=fst_skerry))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color=colors_fst['R'], linetype='F1', size=0.8)+
  geom_point(alpha=0.7, size=2, 
             color=colors_fst['Inversion'], 
             shape=shapes_fst['Inversion']) +
  geom_text(hjust=0.5, vjust=0, aes(label=cp), size=2, face = 'bold', nudge_y = 0.01)+
  annotate(geom="text",x=0.6,y=0.5,
           label=(paste0("R2==",r2_inv)),
           parse=TRUE, color=colors_fst['R'], size=3, fontface = 'bold')+
  ylim(-0.1, 1)+
  xlim(0, 1)+
  xlab('FST in space (Crab-Wave in nearby location)')+
  ylab('FST in time (Skerry 1992 vs 2021)')+
  theme_minimal()+
  #ggtitle(p_title)+
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10),
        #plot.margin = unit(c(0,0.5,0.5,0), "cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.3),
        panel.border = element_rect(colour = "black", fill=NA, size=0.1))


#Generate a PDF file with the correlation plot of Inversions
pdf(width = 5, height = 3, file = './FST_Space_Time/Space_Time_FST_Inversions_2024.pdf')
#Create a grid with two rows: Pie charts and Manhattan plot respectively
plot_grid(plt_inversions, 
          NULL,
          ncol=2,
          nrow=1,
          rel_widths =c(0.73,0.27))
dev.off()
########################################################################################################


########################################################
# Alternatively Plot the two figures one next to the other
########################################################

# extract a legend of the outliers plot which is laid out horizontally
legend_out <- get_legend(
  plt_outliers + theme(legend.position = "bottom")
)

#Hide the legend on outliers plot
plt_outliers = plt_outliers + theme(legend.position="none")

#Create a two-column row with the two plots without legend
prow <- plot_grid(
  plt_outliers,
  plt_inversions,
  nrow = 1,
  ncol = 2
)

#Generate a pdf
pdf(width = 8, height = 4.2, file = './FST_Space_Time/Space_Time_FST_all_2024.pdf')
plot_grid(
  prow,
  legend_out,
  ncol = 1,
  rel_heights = c(0.9, 0.1))
dev.off()
