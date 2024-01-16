################################ FST Trajectories ############################
#Generate the plots that show the trajectories of genetic differentation 
#FST for control (neutal) and spatial outliers 
#in skerry vs Crab and Wave vs Crab over time
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################

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

#Labels to be used in the lagend and title of the figures
name_category_neutral = 'Control'
name_category_outlier = 'Outliers'
name_vs_crab = 'Skerry vs. Crab'
name_vs_wave = 'Skerry vs. Wave'

#Remove SNPs (rows) with missing data in more than 5 individuals in a population
clean_outliers = remove_incomplete_snps(outliers, 5)
clean_neutral = remove_incomplete_snps(neutral, 5)

#Extract the name of the individuals for both control (neutral) and spatial outliers
l_cols_neutral = colnames(clean_neutral)
l_cols_outlier = colnames(clean_outliers)
idx_col_SNP = grep('cat',l_cols) #Index of the column that has the name of the SNP

##############################################
#Subset the data by populations
##############################################

l_cols_crab_1992_neutral = c(idx_col_SNP, grep("DO-92|DO.92", l_cols_neutral))
l_cols_crab_1992_outlier = c(idx_col_SNP, grep("DO-92|DO.92", l_cols_outlier))
crab_1992_neutral = clean_neutral[,l_cols_crab_1992_neutral]
crab_1992_outliers = clean_outliers[,l_cols_crab_1992_outlier]

l_cols_crab_2018_2021_neutral = c(idx_col_SNP, grep("DonCrab|DO-21|DO.21", l_cols_neutral))
l_cols_crab_2018_2021_outlier = c(idx_col_SNP, grep("DonCrab|DO-21|DO.21", l_cols_outlier))
crab_2018_2021_neutral = clean_neutral[,l_cols_crab_2018_2021_neutral]
crab_2018_2021_outliers = clean_outliers[,l_cols_crab_2018_2021_outlier]

l_cols_wave_2018_2021_neutral = c(idx_col_SNP,grep("RefW|W-REF|W.REF", l_cols_neutral))
l_cols_wave_2018_2021_outlier = c(idx_col_SNP,grep("RefW|W-REF|W.REF", l_cols_outlier))
wave_2018_2021_neutral = clean_neutral[,l_cols_wave_2018_2021_neutral]
wave_2018_2021_outliers = clean_outliers[,l_cols_wave_2018_2021_outlier]

l_cols_skerry_2005_neutral = c(idx_col_SNP,grep("SKR-05|SKR.05", l_cols_neutral))
l_cols_skerry_2005_outlier = c(idx_col_SNP,grep("SKR-05|SKR.05", l_cols_outlier))
skerry_2005_neutral = clean_neutral[,l_cols_skerry_2005_neutral]
skerry_2005_outliers = clean_outliers[,l_cols_skerry_2005_outlier]

l_cols_skerry_2018_neutral = c(idx_col_SNP,grep("ExpSkerry", l_cols_neutral))
l_cols_skerry_2018_outlier = c(idx_col_SNP,grep("ExpSkerry", l_cols_outlier))
skerry_2018_neutral = clean_neutral[,l_cols_skerry_2018_neutral]
skerry_2018_outliers = clean_outliers[,l_cols_skerry_2018_outlier]

l_cols_skerry_2021_neutral = c(idx_col_SNP,grep("SRK-21|SRK.21", l_cols_neutral))
l_cols_skerry_2021_outlier = c(idx_col_SNP,grep("SRK-21|SRK.21", l_cols_outlier))
skerry_2021_neutral = clean_neutral[,l_cols_skerry_2021_neutral]
skerry_2021_outliers = clean_outliers[,l_cols_skerry_2021_outlier]


##############################################
# Compute the AF
##############################################
alt_af_crab_1992_neutral = compute_alt_af(crab_1992_neutral)
alt_af_crab_1992_outliers = compute_alt_af(crab_1992_outliers)

alt_af_crab_2018_2021_neutral = compute_alt_af(crab_2018_2021_neutral)
alt_af_crab_2018_2021_outliers = compute_alt_af(crab_2018_2021_outliers)

alt_af_wave_2018_2021_neutral = compute_alt_af(wave_2018_2021_neutral)
alt_af_wave_2018_2021_outliers = compute_alt_af(wave_2018_2021_outliers)

alt_af_skerry_2005_neutral = compute_alt_af(skerry_2005_neutral)
alt_af_skerry_2005_outliers = compute_alt_af(skerry_2005_outliers)

alt_af_skerry_2018_neutral = compute_alt_af(skerry_2018_neutral)
alt_af_skerry_2018_outliers = compute_alt_af(skerry_2018_outliers)

alt_af_skerry_2021_neutral = compute_alt_af(skerry_2021_neutral)
alt_af_skerry_2021_outliers = compute_alt_af(skerry_2021_outliers)


##############################################
# Compute FST
##############################################
#neutral
fst_s1992_crab_neutral = compute_fst(alt_af_crab_1992_neutral, alt_af_crab_2018_2021_neutral)
fst_s2005_crab_neutral = compute_fst(alt_af_skerry_2005_neutral, alt_af_crab_2018_2021_neutral)
fst_s2018_crab_neutral = compute_fst(alt_af_skerry_2018_neutral, alt_af_crab_2018_2021_neutral)
fst_s2021_crab_neutral = compute_fst(alt_af_skerry_2021_neutral, alt_af_crab_2018_2021_neutral)

fst_s1992_wave_neutral = compute_fst(alt_af_crab_1992_neutral, alt_af_wave_2018_2021_neutral)
fst_s2005_wave_neutral = compute_fst(alt_af_skerry_2005_neutral, alt_af_wave_2018_2021_neutral)
fst_s2018_wave_neutral = compute_fst(alt_af_skerry_2018_neutral, alt_af_wave_2018_2021_neutral)
fst_s2021_wave_neutral = compute_fst(alt_af_skerry_2021_neutral, alt_af_wave_2018_2021_neutral)

fst_crab_wave_neutral = compute_fst(alt_af_crab_2018_2021_neutral, alt_af_wave_2018_2021_neutral)


#outliers
fst_s1992_crab_outliers = compute_fst(alt_af_crab_1992_outliers, alt_af_crab_2018_2021_outliers)
fst_s2005_crab_outliers = compute_fst(alt_af_skerry_2005_outliers, alt_af_crab_2018_2021_outliers)
fst_s2018_crab_outliers = compute_fst(alt_af_skerry_2018_outliers, alt_af_crab_2018_2021_outliers)
fst_s2021_crab_outliers = compute_fst(alt_af_skerry_2021_outliers, alt_af_crab_2018_2021_outliers)

fst_s1992_wave_outliers = compute_fst(alt_af_crab_1992_outliers, alt_af_wave_2018_2021_outliers)
fst_s2005_wave_outliers = compute_fst(alt_af_skerry_2005_outliers, alt_af_wave_2018_2021_outliers)
fst_s2018_wave_outliers = compute_fst(alt_af_skerry_2018_outliers, alt_af_wave_2018_2021_outliers)
fst_s2021_wave_outliers = compute_fst(alt_af_skerry_2021_outliers, alt_af_wave_2018_2021_outliers)

fst_crab_wave_outliers = compute_fst(alt_af_crab_2018_2021_outliers, alt_af_wave_2018_2021_outliers)

# Estimate the mean Skerry vs Crab FST and the SD for control loci
tmp_list_crab_neutral = c()
tmp_list_crab_neutral = rbind(tmp_list_crab_neutral, c(1992, mean(fst_s1992_crab_neutral), sd(fst_s1992_crab_neutral), name_category_neutral, name_vs_crab))
tmp_list_crab_neutral = rbind(tmp_list_crab_neutral, c(2005, mean(fst_s2005_crab_neutral), sd(fst_s2005_crab_neutral), name_category_neutral, name_vs_crab))
tmp_list_crab_neutral = rbind(tmp_list_crab_neutral, c(2018, mean(fst_s2018_crab_neutral), sd(fst_s2018_crab_neutral), name_category_neutral, name_vs_crab))
tmp_list_crab_neutral = rbind(tmp_list_crab_neutral, c(2021, mean(fst_s2021_crab_neutral), sd(fst_s2021_crab_neutral), name_category_neutral, name_vs_crab))

# Estimate the mean Skerry vs Crab FST and the SD for outlier loci
tmp_list_crab_outliers = c()
tmp_list_crab_outliers = rbind(tmp_list_crab_outliers, c(1992, mean(fst_s1992_crab_outliers), sd(fst_s1992_crab_outliers), name_category_outlier, name_vs_crab))
tmp_list_crab_outliers = rbind(tmp_list_crab_outliers, c(2005, mean(fst_s2005_crab_outliers), sd(fst_s2005_crab_outliers), name_category_outlier, name_vs_crab))
tmp_list_crab_outliers = rbind(tmp_list_crab_outliers, c(2018, mean(fst_s2018_crab_outliers), sd(fst_s2018_crab_outliers), name_category_outlier, name_vs_crab))
tmp_list_crab_outliers = rbind(tmp_list_crab_outliers, c(2021, mean(fst_s2021_crab_outliers), sd(fst_s2021_crab_outliers), name_category_outlier, name_vs_crab))

# Estimate the mean Skerry vs Wave FST and the SD for control loci
tmp_list_wave_neutral = c()
tmp_list_wave_neutral = rbind(tmp_list_wave_neutral, c(1992, mean(fst_s1992_wave_neutral), sd(fst_s1992_wave_neutral), name_category_neutral, name_vs_wave))
tmp_list_wave_neutral = rbind(tmp_list_wave_neutral, c(2005, mean(fst_s2005_wave_neutral), sd(fst_s2005_wave_neutral), name_category_neutral, name_vs_wave))
tmp_list_wave_neutral = rbind(tmp_list_wave_neutral, c(2018, mean(fst_s2018_wave_neutral), sd(fst_s2018_wave_neutral), name_category_neutral, name_vs_wave))
tmp_list_wave_neutral = rbind(tmp_list_wave_neutral, c(2021, mean(fst_s2021_wave_neutral), sd(fst_s2021_wave_neutral), name_category_neutral, name_vs_wave))

# Estimate the mean Skerry vs Wave FST and the SD for outlier loci
tmp_list_wave_outliers = c()
tmp_list_wave_outliers = rbind(tmp_list_wave_outliers, c(1992, mean(fst_s1992_wave_outliers), sd(fst_s1992_wave_outliers), name_category_outlier, name_vs_wave))
tmp_list_wave_outliers = rbind(tmp_list_wave_outliers, c(2005, mean(fst_s2005_wave_outliers), sd(fst_s2005_wave_outliers), name_category_outlier, name_vs_wave))
tmp_list_wave_outliers = rbind(tmp_list_wave_outliers, c(2018, mean(fst_s2018_wave_outliers), sd(fst_s2018_wave_outliers), name_category_outlier, name_vs_wave))
tmp_list_wave_outliers = rbind(tmp_list_wave_outliers, c(2021, mean(fst_s2021_wave_outliers), sd(fst_s2021_wave_outliers), name_category_outlier, name_vs_wave))
  
#Prepare a df with the Crab vs Wave FST
tmp_list_crab_wave_fst = c()
tmp_list_crab_wave_fst = rbind(tmp_list_crab_wave_fst,c('2018+2021', mean(fst_crab_wave_neutral), sd(fst_crab_wave_neutral), name_category_neutral, 'Crab vs Wave'))
tmp_list_crab_wave_fst = rbind(tmp_list_crab_wave_fst,c('2018+2021', mean(fst_crab_wave_outliers), sd(fst_crab_wave_outliers), name_category_outlier, 'Crab vs Wave'))


#############################################################
# Creates the dataframes with the FST estimates
# and calls the plot function.
# Returns a list with the plots of Skerry and Crab vs Wave
#############################################################
generate_FST_plots <- function(){
  #SKERRY VS CRAB FST PLOT
  df_crab_fst = data.frame(matrix(ncol = 5, nrow = 0))
  df_crab_fst = rbind(df_crab_fst, 
                      tmp_list_crab_neutral,
                      tmp_list_crab_outliers)
  colnames(df_crab_fst) <- c('Year','FST', 'SD', 'Loci', 'Population')
  
  #SKERRY VS WAVE FST PLOT
  df_wave_fst = data.frame(matrix(ncol = 5, nrow = 0))
  df_wave_fst = rbind(df_wave_fst, 
                      tmp_list_wave_neutral,
                      tmp_list_wave_outliers)
  colnames(df_wave_fst) <- c('Year','FST', 'SD', 'Loci', 'Population')
  
  df_skerry_wave_crab_fst = rbind(df_crab_fst, df_wave_fst)
  pl_skerry_crab_wave <- plot_fst(df_skerry_wave_crab_fst, 'FST', FALSE)
  
  #REFERENCE CRAB VS WAVE FST PLOT
  df_ref_crab_wave_fst = as.data.frame(tmp_list_crab_wave_fst)
  colnames(df_ref_crab_wave_fst) <- c('Year','FST', 'SD', 'Loci', 'Population')
  pl_ref_crab_wave <- plot_fst(df_ref_crab_wave_fst, 'C vs W', TRUE)
  
  list_fst = list()
  list_fst[['Skerry']] <- pl_skerry_crab_wave
  list_fst[['Ref']] <- pl_ref_crab_wave
  
  #pdf('FST_Trajectories_combined.pdf', width = 12, height = 3)
  #plot_grid(pl_crab_wave, pl_fake_wave, ncol = 2, nrow = 1, rel_widths = c(0.9,0.1))
  #dev.off()
  
  return(list_fst)
}


######################################################
# Generates a line plot with neutral and outlier trajectories
######################################################
plot_fst <- function(data_fst, p_title, p_hide_labels){
  pl <- ggplot(data_fst, 
         aes(x=if(p_hide_labels){factor(Year)}else{as.numeric(Year)},
             y=as.numeric(FST), 
             colour = factor(Population), 
             group = interaction(Population, Loci),
             shape = factor(Loci),
             linetype = factor(Loci))) + 
    geom_line(linewidth = line_width) + 
    geom_point(size = point_size)+
    theme_minimal()+
    ylim(0,0.23)+
    scale_linetype_manual(values = c(Control = "dashed", Neutral = "dashed", Outliers = "solid")) +
    scale_shape_manual(values = c(Control = 8, Neutral = 8, Outliers = 19)) +
    scale_color_manual(values = c(Wave = "#73BADA", Crab = "#E09157",
                                  'Skerry vs. Wave' = "#73BADA", 'Skerry vs. Crab' = "#E09157",
                                  'W' = '#73BADA',
                                  'Crab vs Wave' = '#7A659E')) +
    theme(legend.position=if(show_legend_geneal){"bottom"}else{'none'}, 
          legend.title=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = 0.3),
          plot.title = element_text(hjust = 0.5, size=12, face='bold'),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10))+
    labs(title = p_title, x="", y = "")
  
  if(p_hide_labels){
    
    pl <- pl+theme(axis.text.y = element_blank())+
      labs(x="", y = element_blank())
  }else{
    pl <- pl+labs(y = "FST")+
      scale_x_continuous( breaks= c(1992,1996,2002,2005,2018,2021) )
  }

  return(pl)
}