################# Spatial Outilers Accumulation in a LG ########################
# Generates a boxplot for each LG to show the expected number of 
# spatial outliers with evidence for selection based on the SNP content.
#
# Figure S17
#
#author: "Diego Garcia"
#date: "2023-05-05"
###############################################################################
  
#Set the working directory
setwd('.')

#Load the required libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)

#Define the path of the input files
snps_selected_outliers_file <- "../Data/SEQSNPTM006_OUTLIER_STATUS_OUT.txt"
swedish_outlier_snp_file <- "../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt"

#Rad files and load datasets
snps_selected_outliers <- read.table(snps_selected_outliers_file, header = T, check.names = FALSE)
outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE) #Outliers with full name and chr info

n_draws <- 1000 #Number of random draws to estimate the expectation

#Subset the status of spatial outliers to keep only those with evidence for selection
true_selelected_outliers <- snps_selected_outliers[snps_selected_outliers$Status, ]

#From the large SNP dataset, keep only the desired columns (contig name, LG, and position in LG)
outliers_LG_data = outliers[,c('contig', 'LG' , 'pos')]
outliers_LG_data$SNP = paste0(outliers_LG_data$contig, '_', outliers_LG_data$pos)

#Discard the SNPs that do not have evidence for selection
true_selelected_outliers_LG_data = outliers_LG_data[outliers_LG_data$SNP %in% true_selelected_outliers$SNP,]


#Estimate the Number of outlier loci outside the expected range (with evidence for selection)
num_draw_elements <- dim(true_selelected_outliers_LG_data)[1]

#Randomly draw (n_draws times) the number of selected outliers 
#from the whole set of spatial outliers
df_random_draws <- draw_randomly(outliers_LG_data, num_draw_elements)

#Summarize the selected outliers and add the info to the draws dataframe
selected_outliers_summary <- true_selelected_outliers_LG_data  %>%
  group_by(LG) %>% 
  summarise(count = n(), type = 'Selected')

#Append the actual count of spatial outliers with evidence for selection to 
#the dataframe of the random draws
df_random_draws = rbind(df_random_draws, selected_outliers_summary)

#Generate a PDF file with the output figure
pdf('SelectedOutliersAccumulation_2024.pdf', width=8, height=4)
# Basic box plot with the random draws
p <- ggplot(df_random_draws[df_random_draws$type == 'Draw',], aes(x=as.factor(LG), y=count), alpha = 0.6) + 
  geom_boxplot(aes(color = type)) +
  scale_color_manual(values = c(Draw="#73BADA", Selected="#E09157"), name='Outliers') +
  theme_minimal() + labs(x = "LG", y = "Count")+
  ylim(0,20)
#Add the actual counts of selected spatial outliers
p <- p + geom_point(data = df_random_draws[df_random_draws$type == 'Selected',], 
                    aes(x=as.factor(LG), y=count, color=type), size = 2)

p

dev.off()


####################################################################
#Outlier loci are expected to be evenly distributed across all LGs.
#Draw multiple times randomly from all outliers, the same proportion of true outliers.
####################################################################
draw_randomly <- function(subject_df, num_draw_elements){
  
  #Create a Empty DataFrame with 0 rows and 3 columns
  df_random_draws = data.frame(matrix(nrow = 0, ncol = 3)) 
  
  #Draw randomlly n times and add the summary to the dataframe
  for(d in seq(1:n_draws)){
    rnd_draw <- sample(seq(1:dim(subject_df)[1]),num_draw_elements,replace=FALSE)
    rnd_draw_data <- subject_df[rnd_draw,]
    rnd_draw_data_summary <- rnd_draw_data  %>%
      group_by(LG) %>% 
      summarise(count = n(), type = 'Draw')
    
    df_random_draws = rbind(df_random_draws, rnd_draw_data_summary)
  }
  
  return(df_random_draws)
  
}