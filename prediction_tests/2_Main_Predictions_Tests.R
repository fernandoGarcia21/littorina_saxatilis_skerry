################# Main Predictions Tests #################################
# Perform different statistical tests to evaluate the accuracy 
# of the predictions at phenotypic level, spatial outliers with evidence 
# for selection, and inversions
#
#author: "Diego Garcia"
#date: "2023-10-10"
##########################################################################

#set the working directory
setwd('.')

#load libraries
library(dplyr)
library(stringr)
source('../generals/GeneralsSkerry.R')

#Define the path and name of the files that contain: 

#phenotypes
phenotypes_file <- "../Data/AllShapePhenotypicData.csv"

#spatial outliers
swedish_outlier_snp_file <- "../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt" 

#The status of the SNPs with respect to the expected range (loci with evidence for selection)
outside_envelope_outliers_file <- "../Data/SEQSNPTM006_OUTLIER_STATUS_OUT.txt"

#Thekaryotypes (clusters in the PCA) of the individuals for each inversion
#This file was generated using the Inversions_Individuals_Karyotypes.R scripts
inversions_karyotypes_file <- "../Data/Individuals_karyotypes_inversions.txt"

#ID of the kariotypes in complex inversions: 
#which IDs are the homokaryotypes and which are the heterokaryotypes
#This file was generated using the Inversions_Individuals_Karyotypes.R scripts
complex_karyotypes_ids_file <- "../Data/Complex_Karyotypes_Ids.txt"

#Define the keywords of the sample names (columns) that will be used in the analysis
sample_names_crab = "DO.92|DonCrab|DO.21"
sample_names_wave = "RefW|W.REF"
sample_names_skerry = "SRK.21"
all_samples = paste(sample_names_crab, sample_names_wave, sample_names_skerry, sep = "|")

#####################################################
#Chi-square test on Genetic data of spatial outliers
#####################################################

#Load the outliers genotypes
swedish_outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE)

#Only keep the SNPid column and the columns of actual genotypes (exclude extra info columns)
l_cols = colnames(swedish_outliers)
l_cols_analysis = grep(all_samples, l_cols)
l_cols_analysis = c( grep('cat',l_cols), l_cols_analysis)
swedish_outliers = swedish_outliers[,l_cols_analysis]

#Change the name of the first column 'cat' for 'Contig', for practicity
l_cols = colnames(swedish_outliers)
colnames(swedish_outliers) <- c('Contig', l_cols[2:length(l_cols)])

#Load the information about the outliers inside/outside the expected range 
outside_envelope_outliers <- read.table(outside_envelope_outliers_file, header = T, check.names = FALSE)
#Only keep the names of the outliers outside the expected range (evidence for selection)
outside_envelope_outliers = outside_envelope_outliers[outside_envelope_outliers$Status == TRUE,]


#Filter the genotypes to keep only outliers outside the expected range
contigs_original = str_split_fixed(swedish_outliers$Contig, '_', 3)[,c(1,2)]
contigs_original = paste(contigs_original[,1],contigs_original[,2], sep = "_")
rows_outliers_outside_envelope = which(contigs_original %in% outside_envelope_outliers$SNP)
new_outliers = swedish_outliers[rows_outliers_outside_envelope, ]


#Add the name of the population (group) to each sample
outliers_t = t(new_outliers[,-1])
colnames(outliers_t) = new_outliers[,1]
ind_names = rownames(outliers_t)
outliers_t = cbind(group=ind_names, outliers_t)
outliers_t[grep(sample_names_crab, ind_names),'group'] = 'Crab'
outliers_t[grep(sample_names_wave, ind_names),'group'] = 'Wave'
outliers_t[grep(sample_names_skerry, ind_names),'group'] = 'Skerry'
outliers_t = as.data.frame(outliers_t)

#Count the number of skerry samples, this value will be added to the summary table
df_n_skerry = data.frame(matrix(ncol = 2, nrow=0))
for (col in 2:dim(outliers_t)[2]){
  tmp_group_snp = outliers_t[,c(1,col)]
  tmp_group_snp = na.omit(tmp_group_snp)
  tmp_table = as.data.frame(table(tmp_group_snp))
  df_n_skerry = rbind(df_n_skerry, data.frame('SNP' = colnames(tmp_table)[2], 
                                              N=sum(tmp_table[tmp_table$group == 'Skerry','Freq'])))
}

#Define the dataframes where the Fisher's exact test results will be arranged
df_fisher_results = data.frame(matrix(nrow = 0, ncol = 10))

#Iterate over all columns (SNPs) in the outliers dataset
for (col in 2:dim(outliers_t)[2]){
  
  #Obtain the genotypes of only one SNP at a time (col)
  tmp_group_snp = outliers_t[,c(1,col)]
  tmp_group_snp = na.omit(tmp_group_snp)
  tmp_snp_name = colnames(tmp_group_snp)[2]
  
  #Subset only the genotypes of the skerry and summarize the counts of each genoytpe
  colnames(tmp_group_snp) = c('group', 'genotype')
  tmp_group_snp_skerry = tmp_group_snp[tmp_group_snp$group == 'Skerry', ]
  tmp_group_snp_crab = tmp_group_snp[tmp_group_snp$group == 'Crab', ]
  tmp_group_snp_wave = tmp_group_snp[tmp_group_snp$group == 'Wave', ]
  
  #Count the number of copies of the alternative and reference alleles (in that order)
  tmp_allele_copies_skerry = count_alleles(tmp_group_snp_skerry)
  tmp_allele_copies_crab = count_alleles(tmp_group_snp_crab)
  tmp_allele_copies_wave = count_alleles(tmp_group_snp_wave)
  
  #Define the contingency tables of observed allele counts between Skerry and Crab, and Skerry and Wave 
  #(first row is the alternative allele values and second the reference values)
  dat_crab <- data.frame(
    "Alt_allele" = c(tmp_allele_copies_skerry[1], tmp_allele_copies_crab[1]),
    "Ref_allele" = c(tmp_allele_copies_skerry[2], tmp_allele_copies_crab[2]),
    row.names = c("Skerry", "Crab"),
    stringsAsFactors = FALSE
  )
  colnames(dat_crab) <- c("Alternative", "Refference")
  
  
  dat_wave <- data.frame(
    "Alt_allele" = c(tmp_allele_copies_skerry[1], tmp_allele_copies_wave[1]),
    "Ref_allele" = c(tmp_allele_copies_skerry[2], tmp_allele_copies_wave[2]),
    row.names = c("Skerry", "Wave"),
    stringsAsFactors = FALSE
  )
  colnames(dat_wave) <- c("Alternative", "Refference")
  
  
  #Perform a Fisher's exact test of independence on the skerry vs Crab and skerry vs Wave
  fisher_sc = fisher.test(dat_crab)
  fisher_sw = fisher.test(dat_wave)
  
  #Combine in a DF the results of Skerry vs Wave and Skerry vs Crab
  df_fisher_results = rbind(df_fisher_results,
                            data.frame('SNP' = tmp_snp_name,
                                       'N' = df_n_skerry[df_n_skerry$SNP == tmp_snp_name, 'N'],
                                       'Skerry_Alt_count' = tmp_allele_copies_skerry[1],
                                       'Skerry_Ref_count' = tmp_allele_copies_skerry[2],
                                       'Wave_Alt_count' = tmp_allele_copies_wave[1],
                                       'Wave_Ref_count' = tmp_allele_copies_wave[2],
                                       'Crab_Alt_count' = tmp_allele_copies_crab[1],
                                       'Crab_Ref_count' = tmp_allele_copies_crab[2],
                                       'P_sw' = format(fisher_sw$p.value, digits = 5),
                                       'P_sc' = format(fisher_sc$p.value, digits = 5)))
  
}

#Sort the dataframe by the p-value of the Skerry vs Wave
df_chi_outliers_merged = df_fisher_results[order(as.numeric(df_fisher_results$P_sw), decreasing = FALSE),]

#Write the Fisher's test to a file
write.table(df_chi_outliers_merged, 'Fisher-outliers-frequency_based.txt', row.names = FALSE, quote = FALSE, sep = '\t')



#####################################################
#Inversions
#####################################################

#Load the karyotypes of the individuals for each inversion
inversions_karyotypes <- read.table(inversions_karyotypes_file, header = T, check.names = FALSE)

#Load the ids of the kariotypes in complex inversions
complex_karyotypes_ids <- read.table(complex_karyotypes_ids_file, header = T, check.names = FALSE)

ind_names = inversions_karyotypes$Individual

use_rows = grep(all_samples, ind_names)

new_karyotypes = inversions_karyotypes[use_rows,]
new_karyotypes$group = NA
new_karyotypes[grep(sample_names_crab, new_karyotypes$Individual),'group'] = 'Crab'
new_karyotypes[grep(sample_names_wave, new_karyotypes$Individual),'group'] = 'Wave'
new_karyotypes[grep(sample_names_skerry, new_karyotypes$Individual),'group'] = 'Skerry'

inversions_names = unique(new_karyotypes$Inversion)
#inversions_names = inversions_names[!inversions_names %in% c("LGC10.2")]

#Define the dataframes where the Fisher's exact test results will be arranged
df_fisher_results_inversions = data.frame(matrix(nrow = 0, ncol = 13))

#Print the number of skerry samples in the inversions
df_n_skerry = data.frame(matrix(ncol = 2, nrow=0))
for (tmp_inv_name in inversions_names){
  tmp_inv_data = new_karyotypes[new_karyotypes$Inversion == tmp_inv_name,c('group','Karyotype')]
  tmp_table = as.data.frame(table(tmp_inv_data))
  df_n_skerry = rbind(df_n_skerry, data.frame('Inversion' = tmp_inv_name, 
                                              N=sum(tmp_table[tmp_table$group == 'Skerry','Freq'])))
}
write.table(df_n_skerry, 'N_estimates_inversions.txt',  quote = FALSE, row.names = FALSE, sep = '\t')

#Perform Fisher's exact test for inversions
for(tmp_inv_name in inversions_names){
  print(tmp_inv_name)
  tmp_inv_data = new_karyotypes[new_karyotypes$Inversion == tmp_inv_name,c('group','Karyotype')]
  
  tmp_group_inv_skerry = tmp_inv_data[tmp_inv_data$group == 'Skerry', ]
  tmp_group_inv_crab = tmp_inv_data[tmp_inv_data$group == 'Crab', ]
  tmp_group_inv_wave = tmp_inv_data[tmp_inv_data$group == 'Wave', ]
  
  #Count the number of copies of the alternative and reference alleles (in that order)
  tmp_arrangement_copies_skerry = count_inversion_arrangements(tmp_inv_name, tmp_group_inv_skerry, complex_karyotypes_ids)
  tmp_arrangement_copies_crab =   count_inversion_arrangements(tmp_inv_name, tmp_group_inv_crab  , complex_karyotypes_ids)
  tmp_arrangement_copies_wave =   count_inversion_arrangements(tmp_inv_name, tmp_group_inv_wave  , complex_karyotypes_ids)
  
  #Define the contingency tables of observed arrangement counts between Skerry and Crab, and Skerry and Wave 
  #(first column is the Arrangement A, second is Arrangement B, and for complex inversions there is arrangement C)
  dat_crab <- data.frame(
    "Arrangement_A" = c(tmp_arrangement_copies_skerry[1], tmp_arrangement_copies_crab[1]),
    "Arrangement_B" = c(tmp_arrangement_copies_skerry[2], tmp_arrangement_copies_crab[2]),
    row.names = c("Skerry", "Crab"),
    stringsAsFactors = FALSE
  )
  colnames(dat_crab) <- c("A", "B")
  
  dat_wave <- data.frame(
    "Arrangement_A" = c(tmp_arrangement_copies_skerry[1], tmp_arrangement_copies_wave[1]),
    "Arrangement_B" = c(tmp_arrangement_copies_skerry[2], tmp_arrangement_copies_wave[2]),
    row.names = c("Skerry", "Wave"),
    stringsAsFactors = FALSE
  )
  colnames(dat_wave) <- c("A", "B")
  
  #Add the a third column with the counts of Kariotype C for complex inversions
  if(tmp_inv_name %in% complex_karyotypes_ids$Inversion){
    dat_crab$C <- c(tmp_arrangement_copies_skerry[3], tmp_arrangement_copies_crab[3])
    dat_wave$C <- c(tmp_arrangement_copies_skerry[3], tmp_arrangement_copies_wave[3])
  }
  
  #Perform a Fisher's exact test of independence on the skerry vs Crab and skerry vs Wave
  fisher_sc = fisher.test(dat_crab)
  fisher_sw = fisher.test(dat_wave)
  
  #Combine in a DF the results of Skerry vs Wave and Skerry vs Crab
  df_fisher_results_inversions = rbind(df_fisher_results_inversions,
                            data.frame('Inversion' = tmp_inv_name,
                                       'N' = df_n_skerry[df_n_skerry$Inversion == tmp_inv_name, 'N'],
                                       'Skerry_A' = tmp_arrangement_copies_skerry[1],
                                       'Skerry_B' = tmp_arrangement_copies_skerry[2],
                                       'Skerry_C' = if(tmp_inv_name %in% complex_karyotypes_ids$Inversion){tmp_arrangement_copies_skerry[3]}else{'-'},
                                       'Wave_A' = tmp_arrangement_copies_wave[1],
                                       'Wave_B' = tmp_arrangement_copies_wave[2],
                                       'Wave_C' = if(tmp_inv_name %in% complex_karyotypes_ids$Inversion){tmp_arrangement_copies_wave[3]}else{'-'},
                                       'Crab_A' = tmp_arrangement_copies_crab[1],
                                       'Crab_B' = tmp_arrangement_copies_crab[2],
                                       'Crab_C' = if(tmp_inv_name %in% complex_karyotypes_ids$Inversion){tmp_arrangement_copies_crab[3]}else{'-'},
                                       'P_sw' = format(fisher_sw$p.value, digits = 5),
                                       'P_sc' = format(fisher_sc$p.value, digits = 5)))
  
}

#Sort the dataframe by the p-value of the Skerry vs Wave
df_fisher_results_inversions = df_fisher_results_inversions[order(as.numeric(df_fisher_results_inversions$P_sw), decreasing = FALSE),]

#Write the Fisher's test to a file
write.table(df_fisher_results_inversions, 'Fisher-inversions-frequency_based.txt', row.names = FALSE, quote = FALSE, sep = '\t')



#####################################################
#Phenotypic data
#####################################################
phenotypes <- read.csv(phenotypes_file, header = T, check.names = FALSE, na.strings="NA", row.names=NULL)

subset_2021 = phenotypes[phenotypes$population %in% c('Skerry 2021',
                                                      'Wave 2018','Wave 2021',
                                                      'Crab 1992','Crab 2018','Crab 2021'),
                         c('population', 'ln(gw)', 'ln(gh)', 
                           'r0', 'h0', 'a0', 'c', 
                           'shell_length', 
                           'ridged', 'colour', 
                           'tesselated', 'avg_thickness')]

subset_2021$group = ''

subset_2021[subset_2021$population == 'Skerry 2021', ]$group = 'Skerry'
subset_2021[subset_2021$population %in% c('Wave 2018', 'Wave 2021'), ]$group = 'Wave'
subset_2021[subset_2021$population %in% c('Crab 1992', 'Crab 2018', 'Crab 2021'), ]$group = 'Crab'

subset_tesselated = subset_2021[subset_2021$tesselated %in% c(0,1), c('group','tesselated')]
subset_ridged = subset_2021[subset_2021$ridged %in% c(0,1), c('group','ridged')]
subset_colour = subset_2021[, c('group','colour')]
#####################################################

##########################################################
#     CHI-SQUARE Test only for qualitative traits
#########################################################
#The Skerry sample size is 44 individuals
#and the Wave ecotype (2021+2018) is 69
#The chi-square test is affected by sample size,
#so we randomly sampled 44 individuals from Wave 100 times
#and estimate the average p values of the 100 t tests
#Set the seed so the random sampling will not change when re-run the code
set.seed(5)
r_replicates = 100
df_x_sqared_sw = data.frame(matrix(nrow = 0, ncol = 3))
df_x_sqared_sc = data.frame(matrix(nrow = 0, ncol = 3))

df_df_sw = data.frame(matrix(nrow = 0, ncol = 3))
df_df_sc = data.frame(matrix(nrow = 0, ncol = 3))

df_pvalues_sw = data.frame(matrix(nrow = 0, ncol = 3))
df_pvalues_sc = data.frame(matrix(nrow = 0, ncol = 3))

for(i in 1:r_replicates){
  
  tmp_chis_tesselated = chi_test(subset_tesselated)
  tmp_chis_ridged = chi_test(subset_ridged)
  tmp_chis_colour = chi_test(subset_colour)
  
  #Merge the results of Skerry vs Wave
  df_x_sqared_sw = rbind(df_x_sqared_sw,
                         data.frame('X2 Ridged' = tmp_chis_ridged[['sw']]$statistic,
                                    'X2 Colour' = tmp_chis_colour[['sw']]$statistic,
                                    'X2 Tesselation'= tmp_chis_tesselated[['sw']]$statistic))
  
  df_df_sw = rbind(df_df_sw,
                   data.frame('DF Ridged' = tmp_chis_ridged[['sw']]$parameter,
                              'DF Colour' = tmp_chis_colour[['sw']]$parameter,
                              'DF Tesselation'= tmp_chis_tesselated[['sw']]$parameter))
  
  df_pvalues_sw = rbind(df_pvalues_sw,
                        data.frame('p-value Ridged' = tmp_chis_ridged[['sw']]$p.value,
                                   'p-value Colour' = tmp_chis_colour[['sw']]$p.value,
                                   'p-value Tesselation'= tmp_chis_tesselated[['sw']]$p.value))
  
  #Merge the results of Skerry vs Crab
  df_x_sqared_sc = rbind(df_x_sqared_sc,
                         data.frame('X2 Ridged' = tmp_chis_ridged[['sc']]$statistic,
                                    'X2 Colour' = tmp_chis_colour[['sc']]$statistic,
                                    'X2 Tesselation'= tmp_chis_tesselated[['sc']]$statistic))
  
  df_df_sc = rbind(df_df_sc,
                   data.frame('DF Ridged' = tmp_chis_ridged[['sc']]$parameter,
                              'DF Colour' = tmp_chis_colour[['sc']]$parameter,
                              'DF Tesselation'= tmp_chis_tesselated[['sc']]$parameter))
  
  df_pvalues_sc = rbind(df_pvalues_sc,
                        data.frame('p-value Ridged' = tmp_chis_ridged[['sc']]$p.value,
                                   'p-value Colour' = tmp_chis_colour[['sc']]$p.value,
                                   'p-value Tesselation'= tmp_chis_tesselated[['sc']]$p.value))
}


print('Chi-square Skerry 2021 vs Wave')
colMeans(df_x_sqared_sw) 
unique(df_df_sw)
colMeans(df_pvalues_sw)

print('Chi-square Skerry 2021 vs Crab')
colMeans(df_x_sqared_sc)
unique(df_df_sc)
colMeans(df_pvalues_sc)



##########################################################
#     T-Test only for quantitative traits
#########################################################

quantitative_traits = c('shell_length', 'avg_thickness', 
                        'ln(gw)', 'ln(gh)', 
                        'r0', 'h0', 'a0', 
                        'c')

df_t_test_sw = data.frame(matrix(nrow = 0, ncol = 4))
df_t_test_sc = data.frame(matrix(nrow = 0, ncol = 4))

for(tmp_qt in quantitative_traits) {
  tmp_subset_qt_Skerry = subset_2021[subset_2021$group == 'Skerry' ,tmp_qt]
  tmp_subset_qt_wave = subset_2021[subset_2021$group == 'Wave' ,tmp_qt]
  tmp_subset_qt_crab = subset_2021[subset_2021$group == 'Crab' ,tmp_qt]
  
  t_test_sw = t.test(tmp_subset_qt_Skerry, tmp_subset_qt_wave)
  t_test_sc = t.test(tmp_subset_qt_Skerry, tmp_subset_qt_crab)
  
  df_t_test_sw = rbind(df_t_test_sw, 
                       data.frame(
                         'Trait' = tmp_qt,
                         'T' = t_test_sw$statistic,
                         'DF' = t_test_sw$parameter,
                         'P' = t_test_sw$p.value
                       ))
  
  df_t_test_sc = rbind(df_t_test_sc, 
                       data.frame(
                         'Trait' = tmp_qt,
                         'T' = t_test_sc$statistic,
                         'DF' = t_test_sc$parameter,
                         'P' = t_test_sc$p.value
                       ))
}
#Show the results of the t-tests on quantitative traits
df_t_test_sw #Skerry vs Wave
df_t_test_sc #Skerry vs Crab



#######################################################
# 
#             FUNCTIONS
#
######################################################

#######################################################
# Estimates the number of copies of the two (simple) 
# or three (complex) arrangements in the sample
# returns a vector with two elements: 
# alt and ref copy number respectively
######################################################
count_inversion_arrangements <- function(p_inv_name, p_individuals_karyotypes, p_complex_karyotype_ids){
  #Generate a summary table with the number of inviduals that have each karyotype in the inversion
  tmp_table = as.data.frame(table(p_individuals_karyotypes))
  
  count_arrangement_A = 0
  count_arrangement_B = 0
  count_arrangement_C = 0
  
  #if the inversion is complex
  if(p_inv_name %in% p_complex_karyotype_ids$Inversion){
    
    #Count the number of copies of the homokaryotypes (A,B,C)
    if(tmp_kariotypes_inv$A %in% tmp_table$Karyotype){
      count_arrangement_A = count_arrangement_A + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$A, 'Freq']*2
    }
    if(tmp_kariotypes_inv$B %in% tmp_table$Karyotype){
      count_arrangement_B = count_arrangement_B + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$B, 'Freq']*2
    }
    if(tmp_kariotypes_inv$C %in% tmp_table$Karyotype){
      count_arrangement_C = count_arrangement_C + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$C, 'Freq']*2
    }
    
    #Count the number of copies of the heterokaryotypes (AC,BC,AB)
    if(tmp_kariotypes_inv$AC %in% tmp_table$Karyotype){
      count_arrangement_A = count_arrangement_A + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$AC, 'Freq']*2
      count_arrangement_C = count_arrangement_C + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$AC, 'Freq']*2
    }
    if(tmp_kariotypes_inv$BC %in% tmp_table$Karyotype){
      count_arrangement_B = count_arrangement_B + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$BC, 'Freq']*2
      count_arrangement_C = count_arrangement_C + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$BC, 'Freq']*2
    }
    if(tmp_kariotypes_inv$AB %in% tmp_table$Karyotype){
      count_arrangement_A = count_arrangement_A + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$AB, 'Freq']*2
      count_arrangement_B = count_arrangement_B + tmp_table[tmp_table$Karyotype == tmp_kariotypes_inv$AB, 'Freq']*2
    }
    
    vec_counts = c(count_arrangement_A, count_arrangement_B, count_arrangement_C)
    
  }else{ #The inversion is simple
    #Karyotypes 1 and 3 are the homokaryotypes in simple inversions
    if(1 %in% tmp_table$Karyotype){
      count_arrangement_A = count_arrangement_A + tmp_table[tmp_table$Karyotype == 1, 'Freq']*2
    }
    
    if(3 %in% tmp_table$Karyotype){
      count_arrangement_B = count_arrangement_B + tmp_table[tmp_table$Karyotype == 3, 'Freq']*2
    }
    #Kariotype 2 is the heterokaryotype in simple inversions
    if(2 %in% tmp_table$Karyotype){
      count_arrangement_A = count_arrangement_A + tmp_table[tmp_table$Karyotype == 2, 'Freq']
      count_arrangement_B = count_arrangement_B + tmp_table[tmp_table$Karyotype == 2, 'Freq']
    }
    
    vec_counts = c(count_arrangement_A, count_arrangement_B)
  }
  
  return( vec_counts )
}

#######################################################
# Estimates the number of copies of the alternative and 
# the refference alleles in the sample
# returns a vector with two elements: 
# alt and ref copy number respectively
######################################################
count_alleles <- function(p_group_snp_genotypes){
  tmp_table = as.data.frame(table(p_group_snp_genotypes))
  
  #Count the number of copies of the alternative and reference alleles in the sample
  count_alt_allele = 0
  count_ref_allele = 0
  
  if('1/1' %in% tmp_table$genotype){
    count_alt_allele = count_alt_allele + tmp_table[tmp_table$genotype == '1/1', 'Freq']*2
  }
  
  if('0/0' %in% tmp_table$genotype){
    count_ref_allele = count_ref_allele + tmp_table[tmp_table$genotype == '0/0', 'Freq']*2
  }
  
  if('0/1' %in% tmp_table$genotype){
    count_alt_allele = count_alt_allele + tmp_table[tmp_table$genotype == '0/1', 'Freq']
    count_ref_allele = count_ref_allele + tmp_table[tmp_table$genotype == '0/1', 'Freq']
  }
  
  if('1/0' %in% tmp_table$genotype){
    count_alt_allele = count_alt_allele + tmp_table[tmp_table$genotype == '1/0', 'Freq']
    count_ref_allele = count_ref_allele + tmp_table[tmp_table$genotype == '1/0', 'Freq']
  }
  
  return( c(count_alt_allele, count_ref_allele) )
}


##########################################
# Sample n elements from Crab and n elements from Wave
# and performns chi-square of Skerry vs Wave and Skerry vs Crab.
# Returns a list with two elements: Crab (sc) and Wave (sw)
##########################################
chi_test <- function(p_subset_2021){
  #Subsample the Wave ecotype to have equal number of individuals as in Skerry
  rows_wave = which(p_subset_2021$group == 'Wave')
  rows_crab = which(p_subset_2021$group == 'Crab')
  n_skerry = length(which(p_subset_2021$group == 'Skerry'))
  out_list = list()
  
  diff_categories_sw = 0
  diff_categories_sc = 0
  
  #The sampling is repeated until we find more than one category in the evaluated variable
  #Chi square cannot test if there is one category 
  aux_cont = 0
  while(diff_categories_sw < 2 & aux_cont < 10){
    subset_wave = p_subset_2021[sample(rows_wave, n_skerry, replace = FALSE), ]
    subset_2021_even_sw = rbind(p_subset_2021[p_subset_2021$group == 'Skerry', ], subset_wave)
    diff_categories_sw = dim(table(subset_2021_even_sw))[2]
    aux_cont = aux_cont+1
  }
  
  if(aux_cont > 1){
    print(paste('SW tried',aux_cont, 'times'))
  }
  
  aux_cont = 0
  while(diff_categories_sc < 2 & aux_cont < 10){
    subset_crab = p_subset_2021[sample(rows_crab, n_skerry, replace = FALSE), ]
    subset_2021_even_sc = rbind(p_subset_2021[p_subset_2021$group == 'Skerry', ], subset_crab)
    diff_categories_sc = dim(table(subset_2021_even_sc))[2]
    aux_cont = aux_cont+1
  }
  
  if(aux_cont > 1){
    print(paste('SC tried',aux_cont, 'times'))
  }
  
  if(diff_categories_sw > 1 & diff_categories_sc > 1){
    CHIS_sw <- chisq.test(subset_2021_even_sw[,1], subset_2021_even_sw[,2])
    CHIS_sc <- chisq.test(subset_2021_even_sc[,1], subset_2021_even_sc[,2])
    
    out_list[['sw']] = CHIS_sw
    out_list[['sc']] = CHIS_sc
  }
  return (out_list)
}
