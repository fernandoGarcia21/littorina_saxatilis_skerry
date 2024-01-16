################# Genomic Distance to FST Correlation ########################
# Generates a histogram with the genome-wide distribution of FST 
# in control loci with respect to the genomic distance to a spatial outlier.
# This plot was used to evaluate signals of hitch-hiking.
#
# Figure S18
#
#author: "Diego Garcia"
#date: "2023-05-22"
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
source('../generals/GeneralsSkerry.R')

#The number of times the genetic map coordinates will be shuffled (yellow lines)
n_shuffles = 100

#Define the path and name of the files that contain control and spatial outliers SNP genotypes
swedish_neutral_snp_file <- "../Data/SkerryExperiment_Neutral_NOLG12_Swedenfilter.txt"
swedish_outlier_snp_file <- "../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt"

#Path to the file with genetic map positions of the contigs
contig_pos_gm_file <- "../Data/contigPositionsGM.txt"

#Read and load the data files with the genotypes
neutral <- read.table(swedish_neutral_snp_file, header = T, check.names = FALSE)
outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE)
pos_gm_data <- read.table(contig_pos_gm_file, header = T, check.names = FALSE)

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

#Change the name of the first column 'cat' for 'Contig', for practicity
l_cols = colnames(neutral)
colnames(neutral) <- c('Contig', l_cols[2:length(l_cols)])
colnames(outliers) <- c('Contig', l_cols[2:length(l_cols)])


#Remove SNPs (rows) with missing data in more than 5 individuals in a population
neutral = remove_incomplete_snps(neutral, 5)
outliers = remove_incomplete_snps(outliers, 5)

#Subset the SNP dataset by skerry population
data_crab_1992 =   get_data_by_population('DO.92',     neutral, outliers, l_cols)
data_skerry_2005 = get_data_by_population("SKR.05",    neutral, outliers, l_cols)
data_skerry_2018 = get_data_by_population("ExpSkerry", neutral, outliers, l_cols)
data_skerry_2021 = get_data_by_population("SRK.21",    neutral, outliers, l_cols)

#Subset the SNP dataset by refference populations
data_crab = get_data_by_population("DonCrab|DO.21",    neutral, outliers, l_cols)
data_wave = get_data_by_population("RefW|W.REF",       neutral, outliers, l_cols)


#Compute AF neutral loci (control)
crab_neutral_af = compute_alt_af(data_crab[[1]])
wave_neutral_af = compute_alt_af(data_wave[[1]])
skerry_neutral_af = list(compute_alt_af(data_crab_1992[[1]]),
                         compute_alt_af(data_skerry_2005[[1]]),
                         compute_alt_af(data_skerry_2018[[1]]),
                         compute_alt_af(data_skerry_2021[[1]]))

#Compute AF spatial outlier loci
crab_outlier_af = compute_alt_af(data_crab[[2]])
wave_outlier_af = compute_alt_af(data_wave[[2]])
skerry_outlier_af = list(compute_alt_af(data_crab_1992[[2]]),
                         compute_alt_af(data_skerry_2005[[2]]),
                         compute_alt_af(data_skerry_2018[[2]]),
                         compute_alt_af(data_skerry_2021[[2]]))

#Compute FST for control loci
fst_s1992_crab_neutral = compute_fst(skerry_neutral_af[[1]], crab_neutral_af)
fst_s2005_crab_neutral = compute_fst(skerry_neutral_af[[2]], crab_neutral_af)
fst_s2018_crab_neutral = compute_fst(skerry_neutral_af[[3]], crab_neutral_af)
fst_s2021_crab_neutral = compute_fst(skerry_neutral_af[[4]], crab_neutral_af)

#Compute FST for spatial outlier loci
fst_s1992_crab_outlier = compute_fst(skerry_outlier_af[[1]], crab_outlier_af)
fst_s2005_crab_outlier = compute_fst(skerry_outlier_af[[2]], crab_outlier_af)
fst_s2018_crab_outlier = compute_fst(skerry_outlier_af[[3]], crab_outlier_af)
fst_s2021_crab_outlier = compute_fst(skerry_outlier_af[[4]], crab_outlier_af)


#Create dataframes with FST estimates for control and spatial outliers
fst_crab_neutral = data.frame('SNP' = neutral$Contig,
                              '1992' = fst_s1992_crab_neutral,
                              '2005' = fst_s2005_crab_neutral,
                              '2018' = fst_s2018_crab_neutral,
                              '2021' = fst_s2021_crab_neutral)
fst_crab_outlier = data.frame('SNP' = outliers$Contig,
                              '1992' = fst_s1992_crab_outlier,
                              '2005' = fst_s2005_crab_outlier,
                              '2018' = fst_s2018_crab_outlier,
                              '2021' = fst_s2021_crab_outlier)


colnames(fst_crab_neutral) = c('SNP','1992','2005','2018','2021')
colnames(fst_crab_outlier) = c('SNP','1992','2005','2018','2021')

# obtain positions in genetic map
df_positions_neutral = generate_SNP_coordinates(neutral)
df_positions_outlier = generate_SNP_coordinates(outliers)



# plot the distribution of FST with respect to genetic distance from the nearest outlier
# Then, shuffle the coordinates and plot the resultant linear model, do it 100 times

#Estimate the distance of a control loci to the nearest outlier
list_dist_from_outl = estimate_distance_from_outlier(df_positions_neutral, df_positions_outlier)
new_positions_neutral = df_positions_neutral[df_positions_neutral$SNP %in% list_dist_from_outl$SNP, ]
new_dist_cM = list_dist_from_outl$dist_cM

#Keep only the FST of control SNPs (neutral) that have info of cM distance to outlier.
#There are neutral loci (e.g. LG 10) that do not have candidate outliers in the LG.
new_fst_neutral = fst_crab_neutral[fst_crab_neutral$SNP %in% list_dist_from_outl$SNP, ]

#Create a PDF file with the output plots
pdf('FST-Distance_shuffled_SupplementaryFigure_Final_2024.pdf', width = 10, height = 10)

par(mfrow = c(2, 2))
plot_distance_stat_shuffles(new_positions_neutral, df_positions_outlier, new_dist_cM, new_fst_neutral[,'1992'], 'cM', 'FST', 'Skerry 1992 vs Crab', '#D35FB7', FALSE, 2)

plot_distance_stat_shuffles(new_positions_neutral, df_positions_outlier, new_dist_cM, new_fst_neutral[,'2005'], 'cM', 'FST', 'Skerry 2005 vs Crab', '#D35FB7', FALSE, 2)

plot_distance_stat_shuffles(new_positions_neutral, df_positions_outlier, new_dist_cM, new_fst_neutral[,'2018'], 'cM', 'FST', 'Skerry 2018 vs Crab', '#D35FB7', FALSE, 2)

plot_distance_stat_shuffles(new_positions_neutral, df_positions_outlier, new_dist_cM, new_fst_neutral[,'2021'], 'cM', 'FST', 'Skerry 2021 vs Crab', '#D35FB7', FALSE, 2)

dev.off()


##############################################################
# Plot FST in Y axis vs binned genetic distance in X axis
#p_cm_bin is the number of centimorgans for each bin
#############################################################
plot_distance_stat_shuffles <- function(df_positions_neutral, 
                                        df_positions_outlier,
                                        dist_x_axis, stat_y_axis,
                                        p_xlab, p_ylab, 
                                        p_title, p_color,
                                        is_scatter,
                                        p_cm_bin){
  
  tmp_positions_neutral_shuffled = df_positions_neutral
  original_neutral_pos = df_positions_neutral$pos
  y_r_label = 0.5
  x_r_label = 30
  
  #Plot either a scatter plot or a binned bar plot
  if(is_scatter){
    plot(dist_x_axis, 
         stat_y_axis, 
         xlab = p_xlab, 
         ylab = p_ylab,
         main = p_title,
         col = alpha('black', 0.3), pch=16)
    
    y_r_label = max(stat_y_axis) * 0.35
    x_r_label = max(dist_x_axis) * 0.85
    
  }else{
    tmp_df_distances = data.frame(dist=dist_x_axis, FST=stat_y_axis)
    #Create barplot
    # With p_cm_bin bins
    # Upper limmit of the last group
    up_lim = ceiling(max(tmp_df_distances$dist))
    up_lim = up_lim + up_lim%%p_cm_bin
    tmp_df_distances$Groups <- cut(x=tmp_df_distances$dist, include.lowest = TRUE, 
                                   breaks=seq(from=0.0, to=up_lim, by = p_cm_bin))
    
    # Estimate mean FST by bins
    summaryFST = tapply(tmp_df_distances$FST, tmp_df_distances$Groups, mean)
    
    #Count the number of SNPs within each bin
    count_snps_bin = tmp_df_distances %>% count(Groups, sort = TRUE)
    df_counts = data.frame(Groups = names(summaryFST), n = 0, FST = as.numeric(summaryFST))
    rownames(df_counts) = df_counts$Groups
    df_counts[as.character(count_snps_bin$Groups), "n"] <- count_snps_bin$n
    df_counts[is.na(df_counts$FST), ]$FST = 0
    
    # Create barplot with average FST for the beans
    bp <- barplot(height = summaryFST, xlab = p_xlab, ylab = p_ylab, main = p_title, cex.axis=0.75, ylim=c(0,0.18), cex.names = 0.5, mgp = c(1, 0.0, 0), col = '#73BADA')
    text(bp, df_counts$FST, df_counts$n ,cex=0.5,pos=3) 
    
    #y_r_label = max(summaryFST[!is.na(summaryFST)]) * 0.5
    y_r_label = 0.08
    x_r_label = 27
  }
  
  #Shuffle the position of neutral loci n times and draw the linear regressions
  set.seed(5) #Seed for the random sampling
  
  for(i in 1:n_shuffles){
    tmp_positions_neutral_shuffled$pos = sample(original_neutral_pos)
    list_dist_from_outl_shuffled = estimate_distance_from_outlier(tmp_positions_neutral_shuffled, df_positions_outlier)
    cM_dist_from_outl_shuffled = list_dist_from_outl_shuffled$dist_cM
    abline(lm(stat_y_axis ~ cM_dist_from_outl_shuffled), col='#FEFE62', lwd = 0.3, alpha=0.6)
  }
  #Draw the linear regression of the actual data
  abline(lm(stat_y_axis ~ dist_x_axis), col=p_color, lwd = 1)
  
  #Estimate the correlation coefficient of x and y values
  cor_v = round(cor(dist_x_axis,stat_y_axis),digits=2)
  text(x=x_r_label, y=y_r_label, paste('r =', cor_v), col=p_color, cex=1)
}


########################################################################
# identify the closest candidate outlier 
# from each neutral locus
########################################################################
estimate_distance_from_outlier <- function(p_df_positions_neutral, p_df_positions_outliers){
  tmp_list_dist_from_outl = data.frame(matrix(nrow = 0,ncol = 2))
  for(i in 1:nrow(p_df_positions_neutral)){
    tmp_neutral = p_df_positions_neutral[i,]
    tmp_chrs_out = p_df_positions_outliers[p_df_positions_outliers$chr == as.numeric(tmp_neutral['chr']),]
    #There are chromosomes with no outliers (e.g. LG 10)
    if(nrow(tmp_chrs_out) > 0){
      tmp_dist = tmp_chrs_out$pos - tmp_neutral$pos
      tmp_nearest_idx = which(abs(tmp_dist) == min(abs(tmp_dist)))[1]
      tmp_nearest_out = tmp_chrs_out[tmp_nearest_idx,]
      
      #Add the distance to the closets outlier
      tmp_dist_from_nearest = abs(tmp_nearest_out$pos - tmp_neutral$pos)
      tmp_list_dist_from_outl = rbind(tmp_list_dist_from_outl, 
                                      data.frame(SNP = c(tmp_neutral$SNP), 
                                                 dist_cM = c(tmp_dist_from_nearest)))
    }
  }
  
  colnames(tmp_list_dist_from_outl) = c('SNP','dist_cM')
  return(tmp_list_dist_from_outl)
}



#############################################
# generate a DF with SNP chr and positions
# dataset_ref is used to extract the names of the contigs 
# and the positions in the genetic map
#############################################
generate_SNP_coordinates <- function(dataset_ref){
  
  tmp_pos_in_map = c()
  tmp_contig_names <- str_split_fixed(dataset_ref$Contig, '_', 4)[,1]
  for(tmp_cont in tmp_contig_names){
    tmp_pos = pos_gm_data[pos_gm_data$contig == tmp_cont,]$av
    tmp_pos_in_map <- append(tmp_pos_in_map, tmp_pos)
  }
  
  ref_names <- as.data.frame(dataset_ref$Contig)
  ref_names <- cbind(ref_names, 
                     as.numeric(str_split_fixed(dataset_ref$Contig, '_', 4)[,3]),
                     tmp_pos_in_map,
                     as.numeric(str_split_fixed(dataset_ref$Contig, '_', 4)[,2]),
                     seq(1,length(dataset_ref$Contig))
  )
  colnames(ref_names) <- c('SNP','chr', 'pos', 'pos_in_contig','row_number')
  
  return(ref_names)
}
