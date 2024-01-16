################# Allele Frequency Shift Toward Wave ##########################
# Visually identify the collinear loci that experienced allele frequency 
# change towards the frequencies in Wave.
#
# Figure S14: Allele frequency scatter plot
# Figure S15: Allele frequency difference AFD scatter plot
#
#author: "Diego Garcia"
#date: "2023-05-05"
###############################################################################

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
source('../generals/GeneralsSkerry.R')

#Define the path and name of the files that contain control and spatial outliers SNP genotypes
swedish_neutral_snp_file <- "../Data/SkerryExperiment_Neutral_NOLG12_Swedenfilter.txt"
swedish_outlier_snp_file <- "../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt"

#Read and load the data files with the genotypes
neutral <- read.table(swedish_neutral_snp_file, header = T, check.names = FALSE)
outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE)
l_cols = colnames(neutral)

#Define the keywords of the sample names (columns) that will be used in the analysis
sample_names_1992 = "DO.92"
sample_names_2005 = "SKR.05"
sample_names_2018 = "DonCrab|ExpSkerry|RefW"
sample_names_2021 = "DO.21|SRK.21|W.REF"

#Merge the column names that will be used in the analysis, this pattern
#is used to subset the columns
sample_names_analysis = paste(sample_names_1992, sample_names_2005, sample_names_2018,sample_names_2021,sep = "|")
#Only keep the SNPid column and the columns of actual genotypes (exclude extra info columns)
l_cols_analysis = grep(sample_names_analysis, l_cols)
l_cols_analysis = c( grep('cat',l_cols), l_cols_analysis)
neutral = neutral[,l_cols_analysis]
outliers = outliers[,l_cols_analysis]

#Remove SNPs (rows) with missing data in more than 5 individuals in a population
neutral = remove_incomplete_snps(neutral, 5)
outliers = remove_incomplete_snps(outliers, 5)
l_cols = colnames(neutral)


#Define the colour and shape palettes of the scatter plots
colors_crab_skerry_wave = c('Wave' = '#73BADA', 'Skerry' = '#7A659E', 'Crab' = '#E09157')
shapes_crab_skerry_wave = c('Wave' = 15, 'Skerry' = 16, 'Crab' = 2)

#Subset the data by popopulation
data_skerry_2005 = get_data_by_population("SKR.05", neutral, outliers, l_cols)
data_skerry_2018 = get_data_by_population("ExpSkerry", neutral, outliers, l_cols)
data_skerry_2021 = get_data_by_population("SRK.21", neutral, outliers, l_cols)
data_crab = get_data_by_population("DO.92|DonCrab|DO.21", neutral, outliers, l_cols)
data_wave = get_data_by_population("RefW|W.REF", neutral, outliers, l_cols)

#Estimate the alternative allele frequency in Crab
crab_neutral_af = compute_alt_af(data_crab[[1]])
crab_outlier_af = compute_alt_af(data_crab[[2]])

#Estimate the alternative allele frequency in Wave
wave_neutral_af = compute_alt_af(data_wave[[1]])
wave_outlier_af = compute_alt_af(data_wave[[2]])

#Estimate the alternative allele frequency in Skerry
skerry_neutral_af_2005 = compute_alt_af(data_skerry_2005[[1]])
skerry_neutral_af_2018 = compute_alt_af(data_skerry_2018[[1]])
skerry_neutral_af_2021 = compute_alt_af(data_skerry_2021[[1]])
skerry_outlier_af_2005 = compute_alt_af(data_skerry_2005[[2]])
skerry_outlier_af_2018 = compute_alt_af(data_skerry_2018[[2]])
skerry_outlier_af_2021 = compute_alt_af(data_skerry_2021[[2]])


#SCATTER PLOT AF SHIFT NEUTRAL
plt_af_skerry_n_2005 = prepare_scatter_plot(crab_neutral_af, skerry_neutral_af_2005, wave_neutral_af, c('','2005',''), '2005 Control', TRUE)
plt_af_skerry_n_2018 = prepare_scatter_plot(crab_neutral_af, skerry_neutral_af_2018, wave_neutral_af, c('','2018',''), '2018 Control', FALSE)
plt_af_skerry_n_2021 = prepare_scatter_plot(crab_neutral_af, skerry_neutral_af_2021, wave_neutral_af, c('','2021',''), '2021 Control', FALSE)

#SCATTER PLOT AF SHIFT OUTLIERS
plt_af_skerry_o_2005 = prepare_scatter_plot(crab_outlier_af, skerry_outlier_af_2005, wave_outlier_af, c('','2005',''), '2005 Spatial Outliers', FALSE)
plt_af_skerry_o_2018 = prepare_scatter_plot(crab_outlier_af, skerry_outlier_af_2018, wave_outlier_af, c('','2018',''), '2018 Spatial Outliers', FALSE)
plt_af_skerry_o_2021 = prepare_scatter_plot(crab_outlier_af, skerry_outlier_af_2021, wave_outlier_af, c('','2021',''), '2021 Spatial Outliers', FALSE)

pdf('WaveAlleleFrequency_Shift_SupplementaryFigure_Final_2024.pdf', width = 10, height = 6)
plot_grid(plt_af_skerry_n_2005, NULL, plt_af_skerry_n_2018, NULL, plt_af_skerry_n_2021,
          plt_af_skerry_o_2005, NULL, plt_af_skerry_o_2018, NULL, plt_af_skerry_o_2021,
          ncol = 5, nrow = 2,
          rel_widths = c(0.298, 0.02, 0.298, 0.02, 0.298))
dev.off()

#SLOPE AFD NEUTRAL
llp_n_2005 <- prepare_scatter_plot_afd(crab_neutral_af, skerry_neutral_af_2005, wave_neutral_af, neutral, 'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2005 Control', TRUE)
llp_n_2018 <- prepare_scatter_plot_afd(crab_neutral_af, skerry_neutral_af_2018, wave_neutral_af, neutral, 'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2018 Control', FALSE)
llp_n_2021 <- prepare_scatter_plot_afd(crab_neutral_af, skerry_neutral_af_2021, wave_neutral_af, neutral, 'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2021 Control', FALSE)

#SLOPE AFD OUTLIERS
llp_o_2005 <- prepare_scatter_plot_afd(crab_outlier_af, skerry_outlier_af_2005, wave_outlier_af, outliers,  'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2005 Spatial Outliers', TRUE)
llp_o_2018 <- prepare_scatter_plot_afd(crab_outlier_af, skerry_outlier_af_2018, wave_outlier_af, outliers,  'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2018 Spatial Outliers', FALSE)
llp_o_2021 <- prepare_scatter_plot_afd(crab_outlier_af, skerry_outlier_af_2021, wave_outlier_af, outliers,  'AFD Wave vs Crab', 'AFD Skerry vs Crab', '2021 Spatial Outliers', FALSE)

pdf('AFD_SkerryVsCrab_SupplementaryFigure_Final_2024.pdf', width = 10, height = 6)
plot_grid(llp_n_2005, NULL, llp_n_2018, NULL, llp_n_2021,
          llp_o_2005, NULL, llp_o_2018, NULL, llp_o_2021, 
          ncol = 5, nrow = 2, 
          rel_widths = c(0.298, 0.02, 0.298, 0.02, 0.298))
dev.off()

######################################################
# Prepare the AF of three populations to be plotted.
# The crab population is the baseline, so, when
# the af of the wave population is smaller than the crab
# flip the AF by subtracting 1. 
# Eg. 1-wave_af; 1-crab_af; 1-skerry_af;
# So, the wave af will be always above the crab baseline
# and we can identify the alleles of skerry that are more
# similar to wave than to crab.
#
# In addition, estimates the AFD of the wave vs crab (X) and skerry vs crab (Y)
######################################################
prepare_scatter_plot_afd <- function(p_crab_af, p_skerry_af, p_wave_af, dataset_ref, pop1, pop2, p_title, p_first){
  
  num_snps = length(p_crab_af)
  count_above_line = 0
  count_below_line = 0
  
  crab_af = p_crab_af[1:length(p_crab_af)]
  skerry_af = p_skerry_af[1:length(p_skerry_af)]
  wave_af = p_wave_af[1:length(p_wave_af)]
  
  #Identify the Wave allele
  for(n in seq(1,num_snps)){
    if(crab_af[n] > wave_af[n]){
      crab_af[n] = 1-crab_af[n];
      wave_af[n] = 1-wave_af[n];
      skerry_af[n] = 1-skerry_af[n];
    }
    wave_af[n] = wave_af[n] - crab_af[n];
    #skerry_af[n] =  -abs(skerry_af[n] - crab_af[n]);
    skerry_af[n] =  skerry_af[n] - crab_af[n];
    
    #Estimate how many dots are plotted above the line (more similar to wave)
    #and how many below the line (more similar to crab)
    if(skerry_af[n] > 0){
      count_above_line = count_above_line + 1
    }else{
      count_below_line = count_below_line + 1
    }
    
  }
  
  #Create the dataframe to plot
  data_afd <- data.frame(
    wave=wave_af, 
    skerry=skerry_af,
    base = rep(0,length(wave_af))
  )
  
  #Create the plot object
  plt <- ggplot(data_afd, aes(x=wave, skerry))+
    geom_point(alpha=0.3, size=1.3, color = 'black', shape = 16) + 
    geom_segment(aes(x=0,xend=1,y=0,yend=0), color=colors_crab_skerry_wave['Crab'], linetype='solid', size=0.8)+
    geom_segment(aes(x=0,xend=1,y=0,yend=1), color=colors_crab_skerry_wave['Wave'], linetype='dotdash', size=0.8)+
    geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color=colors_crab_skerry_wave['Skerry'], linetype='F1', size=0.8)+
    annotate(geom="text",x=0.9,y=0.5,
             label=(paste0("m==",round(coef(lm(data_afd$skerry~data_afd$wave))[2],digits=2))),
             parse=TRUE, color=colors_crab_skerry_wave['Skerry'], size=3, fontface = 'bold')+
    ylim(-0.5, 1)+
    xlim(0.0, 1)+
    theme_minimal()+
    ggtitle(p_title)+
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10),
          plot.margin = unit(c(0,0.5,0.5,0), "cm"),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.1))
  
  if(p_first){
    plt <- plt + 
      xlab(pop1) +
      ylab(pop2)+
      theme(legend.position = c(0.15, 0.85),
            legend.background = element_rect(fill = "white", color = "black", size = 0.1),
            legend.title = element_blank(),
            legend.key.size = unit(2, "mm"),
            axis.title=element_text(size=8))
  }else{
    plt <- plt + 
      xlab(element_blank()) +
      ylab(element_blank())
  }
  
  return(plt)
}


######################################################
# Prepare the AF of three populations to be plotted.
# The crab population is the baseline, so, when
# the af of the wave population is smaller than the crab
# flip the AF by subtracting 1. 
# Eg. 1-wave_af; 1-crab_af; 1-skerry_af;
# So, the wave af will be always above the crab baseline
# and we can identify the alleles of skerry that are more
# similar to wave than to crab.
######################################################
prepare_scatter_plot <- function(p_crab_af, p_skerry_af, p_wave_af, years, p_title, p_first){
  
  num_snps = length(p_crab_af)
  n_skerry_above_crab = 0
  p_skerry_above_crab = 0
  crab_af = p_crab_af[1:length(p_crab_af)]
  skerry_af = p_skerry_af[1:length(p_skerry_af)]
  wave_af = p_wave_af[1:length(p_wave_af)]
  
  for(n in seq(1,num_snps)){
    if(crab_af[n] > wave_af[n]){
      crab_af[n] = 1-crab_af[n];
      wave_af[n] = 1-wave_af[n];
      skerry_af[n] = 1-skerry_af[n];
    }
    
    #Count how many loci have Skerry AF > Crab AF
    if(skerry_af[n] > crab_af[n]){
      n_skerry_above_crab = n_skerry_above_crab + 1
    }
  }
  
  p_skerry_above_crab = (n_skerry_above_crab/num_snps)*100
  plt <- plot_scatter_af_shifts(crab_af, skerry_af, wave_af, crab_af, years, p_skerry_above_crab, p_title, p_first)
  return (plt)
}


######################################################
# Generates a scatter plot with the allele frequencies
# of each of three populations sorted by crab af
######################################################
plot_scatter_af_shifts <- function(crab_data, skerry_data, wave_data, sort_data, years, p_skerry_to_crab, p_title, p_first){
  indexes_sort = sort.int(sort_data, index.return = TRUE)
  sorted_crab = crab_data[indexes_sort$ix]
  sorted_skerry = skerry_data[indexes_sort$ix]
  sorted_wave = wave_data[indexes_sort$ix]
  
  data = data.frame(
    SNP=rep(seq(1,length(sorted_crab)), 3),
    population = c( rep(paste0("Wave",years[3]), length(sorted_wave)),
                    rep("Skerry", length(sorted_skerry)),
                    rep(paste0("Crab",years[1]), length(sorted_crab))
    ),
    AF = c( sorted_wave, 
            sorted_skerry,
            sorted_crab )
  )
  
  plt <- ggplot(data, aes(x=SNP, y=AF, color=population, shape=population)) + 
    theme_minimal() + 
    geom_point(alpha=0.5, size = 1.3) +
    theme(plot.margin = unit(c(0,0.5,0.5,0), "cm"),
          panel.grid.major = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
    scale_color_manual(values=colors_crab_skerry_wave)+
    scale_shape_manual(values=shapes_crab_skerry_wave)+
    annotate("text", x = length(crab_data)*0.85, y = 0.02, size = 3, fontface = 'bold', color = '#4E4D51',
             label = paste("S > C = ",round(p_skerry_to_crab, digits = 2), "%", sep = ""))+
    ggtitle(p_title)+
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 10))+
    xlim(0,length(crab_data)+10)
    
  
  if(p_first){
    plt <- plt + 
      xlab('Locus (Sorted by Crab AF)') +
      ylab('Frequency of the Wave allele')+
      theme(legend.position = c(0.15, 0.85),
            legend.background = element_rect(fill = "white", color = "black", size = 0.1),
            legend.title = element_blank(),
            legend.key.size = unit(2, "mm"),
            axis.title=element_text(size=8))
  }else{
    plt <- plt + 
      xlab(element_blank()) +
      ylab(element_blank()) +
      theme(legend.position='none')
    
  }
  
  
  return(plt)
}


