############## INVERSIONS GENERATE TRAJECTORIES DATA #######################
# Performns a PCA on inversion SNPs, identify the clusters by k-means,
# and estimate the frequency of each karyotype (cluster) within each population.
# Write two output files:
# Inversions_Trajectories_Frequency_Skerry.txt is a table with the arrangement frequencies in Skerry
# Inversions_Trajectories_Frequency_WaveCrab.txt is a table with the arrangement freequencies in Crab and Wave
#
# These output files are used as input for the Figure 3 plot, the inversion frequency trajectories.
# 
# The output files are also the the basis for:
# Table S6: Frequencies of the Wave arrangement in the skerry and reference populations. 
#
# Additionally, this file permits to generate two auxiliar files:
# Inversions_Frequencies_for_Envelope_November2023.txt & Inversions_SampleSize_for_Envelope.txt
# Which are used to estimate the expected range of frequency change for inversions
#
# author: "Diego Garcia"
# date: "12-Dec-2023"
################################################################################


#Set the working directory
setwd('.')

#Load libraries
library("FactoMineR")

#Define the path and name of the files that contain inversion SNP genotypes and coordinates
inv_file <- "../Data/Inv_Genotypes.txt" #SNP genotypes of the inversions
inv_coordinates <- "../Data/Inv_Coordinates.txt" #Coordinates of the inversions (e.g. start and end in the genetic map)
contig_pos_gm_file <- "../Data/contigPositionsGM.txt" #Coordinates of the contigs in the genetic map

#Read and load data
inv_data <- read.table(inv_file, header = T, check.names = FALSE)
inv_coordinates_data <- read.table(inv_coordinates, header = T, check.names = FALSE)
pos_gm_data <- read.table(contig_pos_gm_file, header = T, check.names = FALSE)

l_cols = colnames(inv_data)

#Define the keywords of the sample names (columns) that will be used in the analysis
sample_names_1992 = "DO-92"
sample_names_2005 = "SKR-05"
sample_names_2018 = "DonCrab|ExpSkerry|RefW"
sample_names_2021 = "DO-21|SRK-21|W-REF"

#Merge the column names that will be used in the analysis, this pattern
#is used to subset the columns
sample_names_analysis = paste(sample_names_1992, sample_names_2005, sample_names_2018,sample_names_2021,sep = "|")

#Subset only the columns that will be used in the analysis
l_cols_pca = grep(sample_names_analysis, l_cols)
l_cols_pca = c(1, l_cols_pca)
inv_data_population = inv_data[,l_cols_pca]


# Split the contig names to extract the LG and the coordinates
splt_headers <- as.data.frame(inv_data_population$Contig)
colnames(splt_headers) <- c('Contig')
splt_headers <- cbind(splt_headers, 
                      as.numeric(str_split_fixed(splt_headers$Contig, '_', 4)[,3]),
                      as.numeric(str_split_fixed(splt_headers$Contig, '_', 4)[,2]),
                      pos_gm_data[match(str_split_fixed(splt_headers$Contig, '_', 4)[,1], pos_gm_data$contig), ]$av,
                      seq(1,length(splt_headers$Contig))
)
colnames(splt_headers) <- c('Contig','LG', 'Coordinate', 'pos', 'row_number')


# Extract the SNP names within each inversion within each LG
# LG 1 has two inversions
tmp_pos_lg = inv_coordinates_data[inv_coordinates_data$LG == 1,]
tmp_headers_lg_1_1 = splt_headers[(splt_headers$LG == 1 
                                   & (splt_headers$pos >= tmp_pos_lg[1,]$inner_start  
                                      & splt_headers$pos <= tmp_pos_lg[1,]$inner_end)),]
tmp_headers_lg_1_2 = splt_headers[(splt_headers$LG == 1 
                                   & (splt_headers$pos >= tmp_pos_lg[2,]$inner_start  
                                      & splt_headers$pos <= tmp_pos_lg[2,]$inner_end)),]
tmp_headers_lg_2 = splt_headers[splt_headers$LG == 2,]
tmp_headers_lg_4 = splt_headers[splt_headers$LG == 4,]
tmp_headers_lg_5 = splt_headers[splt_headers$LG == 5,]
tmp_headers_lg_6 = splt_headers[splt_headers$LG == 6,]
#LG 7 has two inversions
tmp_pos_lg = inv_coordinates_data[inv_coordinates_data$LG == 7,]
tmp_headers_lg_7_1 = splt_headers[(splt_headers$LG == 7 
                                   & (splt_headers$pos >= tmp_pos_lg[1,]$inner_start  
                                      & splt_headers$pos <= tmp_pos_lg[1,]$inner_end)),]
tmp_headers_lg_7_2 = splt_headers[(splt_headers$LG == 7 
                                   & (splt_headers$pos >= tmp_pos_lg[2,]$inner_start  
                                      & splt_headers$pos <= tmp_pos_lg[2,]$inner_end)),]
tmp_headers_lg_9 = splt_headers[splt_headers$LG == 9,]
#LG 10 has two inversions
tmp_pos_lg = inv_coordinates_data[inv_coordinates_data$LG == 10,]
tmp_headers_lg_10_1 = splt_headers[(splt_headers$LG == 10 
                                    & (splt_headers$pos >= tmp_pos_lg[1,]$inner_start  
                                       & splt_headers$pos <= tmp_pos_lg[1,]$inner_end)),]
tmp_headers_lg_10_2 = splt_headers[(splt_headers$LG == 10 
                                    & (splt_headers$pos >= tmp_pos_lg[2,]$inner_start  
                                       & splt_headers$pos <= tmp_pos_lg[2,]$inner_end)),]
tmp_headers_lg_11 = splt_headers[splt_headers$LG == 11,]
#LG 14 has two inversions
tmp_pos_lg = inv_coordinates_data[inv_coordinates_data$LG == 14,]
tmp_headers_lg_14_1 = splt_headers[(splt_headers$LG == 14 
                                    & (splt_headers$pos >= tmp_pos_lg[1,]$inner_start  
                                       & splt_headers$pos <= tmp_pos_lg[1,]$inner_end)),]
tmp_headers_lg_14_2 = splt_headers[(splt_headers$LG == 14 
                                    & (splt_headers$pos >= tmp_pos_lg[2,]$inner_start  
                                       & splt_headers$pos <= tmp_pos_lg[2,]$inner_end)),]
tmp_headers_lg_17 = splt_headers[splt_headers$LG == 17,]


# Extract the data within each inversion
tmp_lg_data_1_1  = inv_data_population[tmp_headers_lg_1_1$row_number, ]
tmp_lg_data_1_2  = inv_data_population[tmp_headers_lg_1_2$row_number, ]
tmp_lg_data_2  = inv_data_population[tmp_headers_lg_2$row_number, ]
tmp_lg_data_4  = inv_data_population[tmp_headers_lg_4$row_number, ]
tmp_lg_data_5  = inv_data_population[tmp_headers_lg_5$row_number, ]
tmp_lg_data_6  = inv_data_population[tmp_headers_lg_6$row_number, ]
tmp_lg_data_7_1  = inv_data_population[tmp_headers_lg_7_1$row_number, ]
tmp_lg_data_7_2  = inv_data_population[tmp_headers_lg_7_2$row_number, ]
tmp_lg_data_9  = inv_data_population[tmp_headers_lg_9$row_number, ]
tmp_lg_data_10_1 = inv_data_population[tmp_headers_lg_10_1$row_number, ]
tmp_lg_data_10_2 = inv_data_population[tmp_headers_lg_10_2$row_number, ]
tmp_lg_data_11 = inv_data_population[tmp_headers_lg_11$row_number, ]
tmp_lg_data_14_1 = inv_data_population[tmp_headers_lg_14_1$row_number, ]
tmp_lg_data_14_2 = inv_data_population[tmp_headers_lg_14_2$row_number, ]
tmp_lg_data_17 = inv_data_population[tmp_headers_lg_17$row_number, ]

#Generate a PCA from each inversion
is_replace_na_genotypes = TRUE  #Impute the NA genotypes with the most common genotype in the population
pca_allp_lg_1_1 <- generate_pca_inversion(tmp_lg_data_1_1)
pca_allp_lg_1_2 <- generate_pca_inversion(tmp_lg_data_1_2)
pca_allp_lg_2 <- generate_pca_inversion(tmp_lg_data_2)
pca_allp_lg_4 <- generate_pca_inversion(tmp_lg_data_4)
#pca_allp_lg_5 <- generate_pca_inversion(tmp_lg_data_5)
pca_allp_lg_6 <- generate_pca_inversion(tmp_lg_data_6)
pca_allp_lg_7_1 <- generate_pca_inversion(tmp_lg_data_7_1)
pca_allp_lg_7_2 <- generate_pca_inversion(tmp_lg_data_7_2)
pca_allp_lg_9 <- generate_pca_inversion(tmp_lg_data_9)
pca_allp_lg_10_1 <- generate_pca_inversion(tmp_lg_data_10_1)
pca_allp_lg_10_2 <- generate_pca_inversion(tmp_lg_data_10_2)
pca_allp_lg_11 <- generate_pca_inversion(tmp_lg_data_11)
pca_allp_lg_14_1 <- generate_pca_inversion(tmp_lg_data_14_1)
#pca_allp_lg_14_2 <- generate_pca_inversion(tmp_lg_data_14_2)
pca_allp_lg_17 <- generate_pca_inversion(tmp_lg_data_17)



############ KARYOTYPE FREQUENCY TRAJECTORIES ###############
# Perform the clustering on the result of the PCA and return PC1 and PC2 coordinates
# The second parameter is the number of expected clusters, 
# the third parameter is the number of principal components to estimate the clusters.
# Complex inversions need to estimate the clusters from PC1 and PC2
# the third parameter (merge_skerry_years) indicates how to deal with sample years.
# When True, sampling years are combined by population for Crab, Wave, and Skerry.
# When False, keep the year of the population. e.g. Skerry 2021, Crab 2018
merge_skerry_years = FALSE

coord_2pc_k_1_1  <- generate_kmeans_clustering(pca_allp_lg_1_1,   3, 1, merge_skerry_years)
coord_2pc_k_1_2  <- generate_kmeans_clustering(pca_allp_lg_1_2,   3, 1, merge_skerry_years)
coord_2pc_k_2    <- generate_kmeans_clustering(pca_allp_lg_2,     3, 1, merge_skerry_years)
coord_2pc_k_4    <- generate_kmeans_clustering(pca_allp_lg_4,     3, 1, merge_skerry_years)
coord_2pc_k_6    <- generate_kmeans_clustering(pca_allp_lg_6,     6, 2, merge_skerry_years)
coord_2pc_k_7_1  <- generate_kmeans_clustering(pca_allp_lg_7_1,   3, 1, merge_skerry_years)
coord_2pc_k_7_2  <- generate_kmeans_clustering(pca_allp_lg_7_2,   3, 1, merge_skerry_years)
coord_2pc_k_9    <- generate_kmeans_clustering(pca_allp_lg_9,     3, 1, merge_skerry_years)
coord_2pc_k_10_1 <- generate_kmeans_clustering(pca_allp_lg_10_1,  3, 1, merge_skerry_years)
coord_2pc_k_10_2 <- generate_kmeans_clustering(pca_allp_lg_10_2,  3, 1, merge_skerry_years)
coord_2pc_k_11   <- generate_kmeans_clustering(pca_allp_lg_11,    3, 1, merge_skerry_years)
coord_2pc_k_14_1 <- generate_kmeans_clustering(pca_allp_lg_14_1,  6, 2, merge_skerry_years)
coord_2pc_k_17   <- generate_kmeans_clustering(pca_allp_lg_17,    3, 1, merge_skerry_years)

############ FREQUENCY TRAJECTORIES OF THE ARRANGEMENTS ###############
############### INCLUDE THE NUMBER OF SNPs within each inversion
trajectory_arr_freq_1_1  <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_1_1  ,dim(pca_allp_lg_1_1$var$coord)[1])
trajectory_arr_freq_1_2  <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_1_2  ,dim(pca_allp_lg_1_2$var$coord)[1])
trajectory_arr_freq_2    <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_2    ,dim(pca_allp_lg_2$var$coord)[1])
trajectory_arr_freq_4    <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_4    ,dim(pca_allp_lg_4$var$coord)[1])
trajectory_arr_freq_6    <- generate_complex_inversion_arrangement_tajectory(coord_2pc_k_6   ,dim(pca_allp_lg_6$var$coord)[1])
trajectory_arr_freq_7_1  <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_7_1  ,dim(pca_allp_lg_7_1$var$coord)[1])
trajectory_arr_freq_7_2  <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_7_2  ,dim(pca_allp_lg_7_2$var$coord)[1])
trajectory_arr_freq_9    <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_9    ,dim(pca_allp_lg_9$var$coord)[1])
trajectory_arr_freq_10_1 <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_10_1 ,dim(pca_allp_lg_10_1$var$coord)[1])
trajectory_arr_freq_10_2 <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_10_2 ,dim(pca_allp_lg_10_2$var$coord)[1])
trajectory_arr_freq_11   <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_11   ,dim(pca_allp_lg_11$var$coord)[1])
trajectory_arr_freq_14_1 <- generate_complex_inversion_arrangement_tajectory(coord_2pc_k_14_1,dim(pca_allp_lg_14_1$var$coord)[1])
trajectory_arr_freq_17   <- generate_simple_inversion_arrangement_tajectory(coord_2pc_k_17   ,dim(pca_allp_lg_17$var$coord)[1])

########################## WRITE TRAJECTORIES INFO FOR DRIFT PLOT ###########

#Add the inversion name to the output dataset
trajectory_arr_freq_1_1  [['Skerry']] $Inversion   <- 'LGC1.1'   
trajectory_arr_freq_1_2  [['Skerry']] $Inversion   <- 'LGC1.2'   
trajectory_arr_freq_2    [['Skerry']]   $Inversion <- 'LGC2.1'   
trajectory_arr_freq_4    [['Skerry']]   $Inversion <- 'LGC4.1'   
trajectory_arr_freq_6    [['Skerry']]   $Inversion <- 'LGC6.1/2' 
trajectory_arr_freq_7_1  [['Skerry']] $Inversion   <- 'LGC7.1'   
trajectory_arr_freq_7_2  [['Skerry']] $Inversion   <- 'LGC7.2'   
trajectory_arr_freq_9    [['Skerry']]   $Inversion <- 'LGC9.1'   
trajectory_arr_freq_10_1 [['Skerry']]$Inversion    <- 'LGC10.1'  
trajectory_arr_freq_10_2 [['Skerry']]$Inversion    <- 'LGC10.2'  
trajectory_arr_freq_11   [['Skerry']]  $Inversion  <- 'LGC11.1'  
trajectory_arr_freq_14_1 [['Skerry']]$Inversion    <- 'LGC14.1/2'
trajectory_arr_freq_17   [['Skerry']]  $Inversion  <- 'LGC17.1'  

trajectory_arr_freq_1_1  [['WaveCrab']] $Inversion   <- 'LGC1.1'   
trajectory_arr_freq_1_2  [['WaveCrab']] $Inversion   <- 'LGC1.2'   
trajectory_arr_freq_2    [['WaveCrab']]   $Inversion <- 'LGC2.1'   
trajectory_arr_freq_4    [['WaveCrab']]   $Inversion <- 'LGC4.1'   
trajectory_arr_freq_6    [['WaveCrab']]   $Inversion <- 'LGC6.1/2' 
trajectory_arr_freq_7_1  [['WaveCrab']] $Inversion   <- 'LGC7.1'   
trajectory_arr_freq_7_2  [['WaveCrab']] $Inversion   <- 'LGC7.2'   
trajectory_arr_freq_9    [['WaveCrab']]   $Inversion <- 'LGC9.1'   
trajectory_arr_freq_10_1 [['WaveCrab']]$Inversion    <- 'LGC10.1'  
trajectory_arr_freq_10_2 [['WaveCrab']]$Inversion    <- 'LGC10.2'  
trajectory_arr_freq_11   [['WaveCrab']]  $Inversion  <- 'LGC11.1'  
trajectory_arr_freq_14_1 [['WaveCrab']]$Inversion    <- 'LGC14.1/2'
trajectory_arr_freq_17   [['WaveCrab']]  $Inversion  <- 'LGC17.1'  

#### Write a file with the trajectories of the frequencies in Skerry
merge_freq_inversions <- rbind(trajectory_arr_freq_1_1   [['Skerry']], 
                               trajectory_arr_freq_1_2 [['Skerry']], 
                               trajectory_arr_freq_2   [['Skerry']],   
                               trajectory_arr_freq_4   [['Skerry']],  
                               trajectory_arr_freq_6   [['Skerry']][,c(1,2,4,5,6,7)],
                               trajectory_arr_freq_7_1 [['Skerry']], 
                               trajectory_arr_freq_7_2 [['Skerry']], 
                               trajectory_arr_freq_9   [['Skerry']],  
                               trajectory_arr_freq_10_1[['Skerry']],
                               trajectory_arr_freq_10_2[['Skerry']],
                               trajectory_arr_freq_11  [['Skerry']],  
                               trajectory_arr_freq_14_1[['Skerry']][,c(1,2,4,5,6,7)],
                               trajectory_arr_freq_17  [['Skerry']])

#This output file is also the the basis for the Skerry columns in 
# Table S6: Frequencies of the Wave arrangement in the skerry and reference populations. 
write.table(merge_freq_inversions, '../Data/Inversions_Trajectories_Frequency_Skerry.txt', append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#### Write a file with the frequencies of the arrangements in the Wave and Crab ecotype
merge_freq_inversions <- rbind(trajectory_arr_freq_1_1  [['WaveCrab']], 
                               trajectory_arr_freq_1_2  [['WaveCrab']], 
                               trajectory_arr_freq_2    [['WaveCrab']],   
                               trajectory_arr_freq_4    [['WaveCrab']],  
                               trajectory_arr_freq_6    [['WaveCrab']][,c(1,2,4,5,6,7)],
                               trajectory_arr_freq_7_1  [['WaveCrab']], 
                               trajectory_arr_freq_7_2  [['WaveCrab']], 
                               trajectory_arr_freq_9    [['WaveCrab']],  
                               trajectory_arr_freq_10_1 [['WaveCrab']],
                               trajectory_arr_freq_10_2 [['WaveCrab']],
                               trajectory_arr_freq_11   [['WaveCrab']],  
                               trajectory_arr_freq_14_1 [['WaveCrab']][,c(1,2,4,5,6,7)],
                               trajectory_arr_freq_17   [['WaveCrab']])

#This output file is also the the basis for the pW2018+2021	& pC2018+2021 columns in 
# Table S6: Frequencies of the Wave arrangement in the skerry and reference populations.
write.table(merge_freq_inversions, '../Data/Inversions_Trajectories_Frequency_WaveCrab.txt', append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)





#################### ADDITIONAL FILES WHICH CONTENT WE USE TO ESTIMATE THE EXPECTED RANGE 

###### Write a file with the frequencies of the inversions in the format
# that is required for neutral envelope analysis
df_invp_envelope = rbind(arrange_inv_p_for_envelope(trajectory_arr_freq_1_1  , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_1_2  , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_2    , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_4    , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_6    , TRUE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_7_1  , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_7_2  , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_9    , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_10_1 , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_10_2 , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_11   , FALSE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_14_1 , TRUE),
                         arrange_inv_p_for_envelope(trajectory_arr_freq_17   , FALSE)
)

write.table(df_invp_envelope, '../Data/Inversions_Frequencies_for_Envelope_November2023.txt', append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)



#Generate a file with the sample size for each inversion, and for each arrangement
df_sample_size_inversions = data.frame(matrix(nrow=0,ncol=3))
df_sample_size_inversions = rbind(df_sample_size_inversions, 
                                  coord_2pc_k_1_1  %>% count(Inversion='LGC1.1'   , Population, sort = TRUE),
                                  coord_2pc_k_1_2  %>% count(Inversion='LGC1.2'   , Population, sort = TRUE),
                                  coord_2pc_k_2    %>% count(Inversion='LGC2.1'   , Population, sort = TRUE),
                                  coord_2pc_k_4    %>% count(Inversion='LGC4.1'   , Population, sort = TRUE),
                                  coord_2pc_k_6    %>% count(Inversion='LGC6.1/2' , Population, sort = TRUE),
                                  coord_2pc_k_7_1  %>% count(Inversion='LGC7.1'   , Population, sort = TRUE),
                                  coord_2pc_k_7_2  %>% count(Inversion='LGC7.2'   , Population, sort = TRUE),
                                  coord_2pc_k_9    %>% count(Inversion='LGC9.1'   , Population, sort = TRUE),
                                  coord_2pc_k_10_1 %>% count(Inversion='LGC10.1'  , Population, sort = TRUE),
                                  coord_2pc_k_10_2 %>% count(Inversion='LGC10.2'  , Population, sort = TRUE),
                                  coord_2pc_k_11   %>% count(Inversion='LGC11.1'  , Population, sort = TRUE),
                                  coord_2pc_k_14_1 %>% count(Inversion='LGC14.1/2', Population, sort = TRUE),
                                  coord_2pc_k_17   %>% count(Inversion='LGC17.1'  , Population, sort = TRUE))

write.table(df_sample_size_inversions, '../Data/Inversions_SampleSize_for_Envelope.txt', append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
########################################################################################################


################################################################################
# Re arrange the structure of the data of each inversion to the
# structure (Inv_ID		pC1992		pW0000	pS2005		pS2018		pS2021).
# that is used to estimate the expected range of frequency change of inversions
################################################################################
arrange_inv_p_for_envelope <- function(p_df_invesion, p_is_complex){
  if(p_is_complex){
    tmp_pskerry = p_df_invesion   [['Skerry']]
    tmp_transpose_freq = rbind(tmp_pskerry[tmp_pskerry$Allele == 'Top 1', ]$Frequency,
                               tmp_pskerry[tmp_pskerry$Allele == 'Top 2', ]$Frequency)
  }else{
    tmp_transpose_freq = t(p_df_invesion   [['Skerry']]['Frequency'])
  }
  
  tmp_transpose_freq = cbind(unique(p_df_invesion   [['Skerry']]['Inversion']), tmp_transpose_freq)
  tmp_pwave = p_df_invesion[['WaveCrab']]
  tmp_pwave = tmp_pwave[tmp_pwave$Year == 'Wave', ]$Frequency
  tmp_all_p_inversion = cbind(tmp_transpose_freq[,1:2], tmp_pwave, tmp_transpose_freq[,3:5])
  colnames(tmp_all_p_inversion) <- c('Inv_ID','pC1992','pW0000','pS2005','pS2018','pS2021')
  
  return(tmp_all_p_inversion)
}




########################################################################################
# ESTIMATES THE FREQUENCY OF THE ARRANGEMENTS WITHIN A COMPLEX INVERSION. 
# COMPLEX INVERSIONS HAVE 3 ARRANGEMTS (A|B|C) AND 3 HETEROKARIOTYPES (AC / BC / AB)
# This function returns a dataframe with the frequency of the arrangments A, B and C
#######################################################################################
estimate_complex_inversion_arrangement_frequencies <- function(p_haplotype_data, p_complex_karyotypes){
  
  #Create an empty dataframe were we will store the info of the frequencies
  inversion_allele_frequencies <- data.frame(Population=character(),
                                             arrangement=integer(), 
                                             count=integer(), 
                                             frequency=integer())
  
  total_copies = sum(p_haplotype_data$count)*2
  pop_name = max(p_haplotype_data$Population)
  
  #When the inversion is complex, we estimate the frequency
  #of each of three arrangements, starting with the homokaryotypes (x2 ALLES).
  count_alleles_AA <- p_haplotype_data[p_haplotype_data$NewGroup == p_complex_karyotypes["A"],]$count * 2
  count_alleles_BB <- p_haplotype_data[p_haplotype_data$NewGroup == p_complex_karyotypes["B"],]$count * 2
  count_alleles_CC <- p_haplotype_data[p_haplotype_data$NewGroup == p_complex_karyotypes["C"],]$count * 2
  
  #The heterokaryotypes only have one copy of each allele
  count_alleles_AC <- p_haplotype_data[p_haplotype_data$NewGroup  == p_complex_karyotypes["AC"],]$count
  count_alleles_BC <- p_haplotype_data[p_haplotype_data$NewGroup  == p_complex_karyotypes["BC"],]$count
  count_alleles_AB <- p_haplotype_data[p_haplotype_data$NewGroup  == p_complex_karyotypes["AB"],]$count
  
  count_A = (if(length(count_alleles_AA) > 0) count_alleles_AA else 0) + 
    (if(length(count_alleles_AC) > 0) count_alleles_AC else 0) + 
    (if(length(count_alleles_AB) > 0) count_alleles_AB else 0)
  
  count_B = (if(length(count_alleles_BB) > 0) count_alleles_BB else 0) + 
    (if(length(count_alleles_BC) > 0) count_alleles_BC else 0) + 
    (if(length(count_alleles_AB) > 0) count_alleles_AB else 0)
  
  count_C = (if(length(count_alleles_CC) > 0) count_alleles_CC else 0) + 
    (if(length(count_alleles_AC) > 0) count_alleles_AC else 0) + 
    (if(length(count_alleles_BC) > 0) count_alleles_BC else 0)
  
  inversion_allele_frequencies = rbind(inversion_allele_frequencies, 
                                       data.frame(Population = pop_name, 
                                                  arrangement = unname(unlist(p_complex_karyotypes["A"])), 
                                                  count = count_A, 
                                                  frequency = count_A / total_copies))
  
  inversion_allele_frequencies = rbind(inversion_allele_frequencies, 
                                       data.frame(Population = pop_name, 
                                                  arrangement = unname(unlist(p_complex_karyotypes["B"])), 
                                                  count = count_B, 
                                                  frequency = count_B / total_copies))
  
  inversion_allele_frequencies = rbind(inversion_allele_frequencies, 
                                       data.frame(Population = pop_name, 
                                                  arrangement = unname(unlist(p_complex_karyotypes["C"])), 
                                                  count = count_C, 
                                                  frequency = count_C / total_copies))
  
  
  return (inversion_allele_frequencies)
}




############################################
# ESTIMATES THE FREQUENCY OF THE ARRANGEMENTS
# WITHIN AN INVERSION. 
# SIMPLE INVERSIONS HAVE 2 ARRANGEMTS (1|0) AND 3 GENOTYPES (11 / 10 / 00)
###########################################
estimate_simple_inversion_arrangement_frequencies <- function(p_haplotype_data){
  
  #Create an empty dataframe were we will store the info of the frequencies
  inversion_allele_frequencies <- data.frame(Population=character(),
                                             arrangement=integer(), 
                                             count=integer(), 
                                             frequency=integer()) 
  
  total_copies = sum(p_haplotype_data$count)*2
  pop_name = max(p_haplotype_data$Population)
  
  #When the inversion is simple, we estimate the frequency
  #of the each of the two possible arrangements (1 or 0) = homokaryotypes 1 and 3
  count_alleles_heterokaryotype = p_haplotype_data[p_haplotype_data$NewGroup == 2,]$count
  if(length(count_alleles_heterokaryotype) == 0){
    count_alleles_heterokaryotype = 0
  }
  for (g in c(1,3)) {
    tmp_karyotype = p_haplotype_data[p_haplotype_data$NewGroup == g,]
    tmp_arrangement_count = 0
    #If there is information for the karyotype g
    if(dim(tmp_karyotype)[1] > 0){
      tmp_karyotype_count = tmp_karyotype$count
      tmp_arrangement_count = (tmp_karyotype_count * 2) + count_alleles_heterokaryotype
    }else{
      tmp_arrangement_count = count_alleles_heterokaryotype
    }
    tmp_arrangement_frequency = tmp_arrangement_count / total_copies
    new_arrangement_info = data.frame(Population = pop_name, 
                                      arrangement = g, 
                                      count = tmp_arrangement_count, 
                                      frequency = tmp_arrangement_frequency)
    inversion_allele_frequencies = rbind(inversion_allele_frequencies, new_arrangement_info)
  }
  
  return (inversion_allele_frequencies)
}


#########################################################################
# CREATES A DATAFRAME WITH THE FREQUENCY OF THE TOP 2 ARRANGEMENTS
# IN THE WAVE ECOTYPE, IN THE DIFFERENT YEARS OF THE SKERRY INCLUDING CRAB 1992
#########################################################################
generate_complex_inversion_arrangement_tajectory <- function(pca_coord_2pc, num_snps_in_inversion){
  aux_name_top_1 = 'Top 1'
  aux_name_top_2 = 'Top 2'
  
  #Create an empty dataframe for the frequency of the two most common arangements in wave
  df_arrangements_frequencies <- data.frame(Year=character(),
                                            Frequency=double(), 
                                            Allele=character(), 
                                            Population=character(),
                                            Color = character())
  
  #Group the pca data of the inversion by population and cluster
  #to estimate the number of snails within each cluster by population
  frequency_df <- pca_coord_2pc %>% 
    group_by(Population, NewGroup) %>% 
    summarise(count = n())
  
  #Identify what clusters are the homokaryotypes and what are the heterokaryotypes
  list_complex_karyotypes = determine_complex_karyotypes(pca_coord_2pc)
  aux_homokaryotypes = unname(unlist(list_complex_karyotypes[c("A","B","C")]))
  
  # Identify the first two arrangements with the largest frequency in wave.
  # Sometimes there is only one arrangement in the population,
  # this happens if there are not homokaryotypes.
  tmp_wave_data <- frequency_df[frequency_df$Population == 'Wave',]
  tmp_wave_arrf_data   <- estimate_complex_inversion_arrangement_frequencies(tmp_wave_data, list_complex_karyotypes)
  tmp_wave_arrf_data <- tmp_wave_arrf_data[order(tmp_wave_arrf_data$count, decreasing = TRUE),]
  
  #Etimate the frequency of the homokariotypes in Crab
  tmp_crab_data <- frequency_df[frequency_df$Population == 'Crab',]
  tmp_crab_arrf_data   <- estimate_complex_inversion_arrangement_frequencies(tmp_crab_data, list_complex_karyotypes)
  tmp_crab_arrf_data <- tmp_crab_arrf_data[order(tmp_crab_arrf_data$count, decreasing = TRUE),]
  
  #Identify the most common arrangement in wave than in crab.
  if(tmp_wave_arrf_data[1,'frequency'] > tmp_crab_arrf_data[1,'frequency']){
    tmp_max_homok_wave = tmp_wave_arrf_data[c(1),]
  }else{
    if(tmp_wave_arrf_data[2,'frequency'] > tmp_crab_arrf_data[2,'frequency']){
      tmp_max_homok_wave = tmp_wave_arrf_data[c(2),]
    }else{
      tmp_max_homok_wave = tmp_wave_arrf_data[c(3),]
    }
  }
  #Add a second row with the allele that is more common in crab than in wave 
  # except the arrangement that was already selected above
  tmp_crab_arrf_data <- tmp_crab_arrf_data[tmp_crab_arrf_data$arrangement != tmp_max_homok_wave$arrangement,]
  tmp_max_homok_wave = rbind(tmp_max_homok_wave, tmp_crab_arrf_data[1,])
  
  #Extract the data of the different years of the skerry population
  tmp_crab_1992_data <- frequency_df[frequency_df$Population == 'Crab 1992',]
  tmp_skerry_2005_data <- frequency_df[frequency_df$Population == 'Skerry 2005',]
  tmp_skerry_2018_data <- frequency_df[frequency_df$Population == 'Skerry 2018',]
  tmp_skerry_2021_data <- frequency_df[frequency_df$Population == 'Skerry 2021',]
  
  #Estimate the frequency of the crab and wave inversion in skerry for different years
  tmp_crab_1992_arrf_data   <- estimate_complex_inversion_arrangement_frequencies(tmp_crab_1992_data, list_complex_karyotypes)
  tmp_skerry_2005_arrf_data <- estimate_complex_inversion_arrangement_frequencies(tmp_skerry_2005_data, list_complex_karyotypes)
  tmp_skerry_2018_arrf_data <- estimate_complex_inversion_arrangement_frequencies(tmp_skerry_2018_data, list_complex_karyotypes)
  tmp_skerry_2021_arrf_data <- estimate_complex_inversion_arrangement_frequencies(tmp_skerry_2021_data, list_complex_karyotypes)
  
  #Extract only the frequencies of the top 2 wave alleles
  tmp_f_1992_wave_alle = tmp_crab_1992_arrf_data  [tmp_crab_1992_arrf_data$arrangement  %in% tmp_max_homok_wave$arrangement, c("arrangement","frequency")]
  tmp_f_2005_wave_alle = tmp_skerry_2005_arrf_data[tmp_skerry_2005_arrf_data$arrangement  %in% tmp_max_homok_wave$arrangement, c("arrangement","frequency")]
  tmp_f_2018_wave_alle = tmp_skerry_2018_arrf_data[tmp_skerry_2018_arrf_data$arrangement  %in% tmp_max_homok_wave$arrangement, c("arrangement","frequency")]
  tmp_f_2021_wave_alle = tmp_skerry_2021_arrf_data[tmp_skerry_2021_arrf_data$arrangement  %in% tmp_max_homok_wave$arrangement, c("arrangement","frequency")]
  
  #Sort the the frequencies according to the order of the top 2 arrangements in wave
  tmp_f_1992_wave_alle = tmp_f_1992_wave_alle[order(match(tmp_f_1992_wave_alle$arrangement,tmp_max_homok_wave$arrangement)),"frequency"]
  tmp_f_2005_wave_alle = tmp_f_2005_wave_alle[order(match(tmp_f_2005_wave_alle$arrangement,tmp_max_homok_wave$arrangement)),"frequency"]
  tmp_f_2018_wave_alle = tmp_f_2018_wave_alle[order(match(tmp_f_2018_wave_alle$arrangement,tmp_max_homok_wave$arrangement)),"frequency"]
  tmp_f_2021_wave_alle = tmp_f_2021_wave_alle[order(match(tmp_f_2021_wave_alle$arrangement,tmp_max_homok_wave$arrangement)),"frequency"]
  
  #Add the top 1 of the wave alleles to the dataframe
  tmp_row_1992_wave_allele_1 = data.frame(Year = 'C 1992', Frequency = tmp_f_1992_wave_alle[1], Allele = aux_name_top_1, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[1])
  tmp_row_2005_wave_allele_1 = data.frame(Year = 'S 2005', Frequency = tmp_f_2005_wave_alle[1], Allele = aux_name_top_1, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[1])
  tmp_row_2018_wave_allele_1 = data.frame(Year = 'S 2018', Frequency = tmp_f_2018_wave_alle[1], Allele = aux_name_top_1, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[1])
  tmp_row_2021_wave_allele_1 = data.frame(Year = 'S 2021', Frequency = tmp_f_2021_wave_alle[1], Allele = aux_name_top_1, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[1])
  
  df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies,
                                             tmp_row_1992_wave_allele_1,
                                             tmp_row_2005_wave_allele_1,
                                             tmp_row_2018_wave_allele_1,
                                             tmp_row_2021_wave_allele_1)
  
  #Add the top 2 of the wave alleles to the dataframe if this the top 2 exists
  if(!is.na(tmp_f_1992_wave_alle[2])){
    tmp_row_aux_wave_allele_2 = data.frame(Year = 'C 1992', Frequency = tmp_f_1992_wave_alle[2], Allele = aux_name_top_2, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[2])
    df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies_skerry,
                                               tmp_row_aux_wave_allele_2)
  }
  
  if(!is.na(tmp_f_2005_wave_alle[2])){
    tmp_row_aux_wave_allele_2 = data.frame(Year = 'S 2005',Frequency = tmp_f_2005_wave_alle[2], Allele = aux_name_top_2, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[2])
    df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies_skerry,
                                               tmp_row_aux_wave_allele_2)
  }
  
  if(!is.na(tmp_f_2018_wave_alle[2])){
    tmp_row_aux_wave_allele_2 = data.frame(Year = 'S 2018', Frequency = tmp_f_2018_wave_alle[2], Allele = aux_name_top_2, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[2])
    df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies_skerry,
                                               tmp_row_aux_wave_allele_2)
  }
  
  if(!is.na(tmp_f_2021_wave_alle[2])){
    tmp_row_aux_wave_allele_2 = data.frame(Year = 'S 2021', Frequency = tmp_f_2021_wave_alle[2], Allele = aux_name_top_2, Population = 'Wave', Color = tmp_max_homok_wave$arrangement[2])
    df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies_skerry,
                                               tmp_row_aux_wave_allele_2)
  }
  
  
  #Obtain the frequency information of the top two arrangements, in the reference populations
  aux_wave_freqs_df = tmp_wave_arrf_data[tmp_wave_arrf_data$arrangement %in% tmp_max_homok_wave$arrangement,]
  aux_crab_freqs_df = tmp_crab_arrf_data[tmp_crab_arrf_data$arrangement %in% tmp_max_homok_wave$arrangement,]
  
  #Check if crab is missing one arrangement present in wave
  aux_arranget_missing = aux_wave_freqs_df[!aux_wave_freqs_df$arrangement %in% aux_crab_freqs_df$arrangement, 'arrangement']
  if(length(aux_arranget_missing) > 0){
    aux_crab_freqs_df = rbind(aux_crab_freqs_df,
                              data.frame(Population = 'Crab',
                                         arrangement = aux_arranget_missing,
                                         count = 0,
                                         frequency = 0))
  }
  
  #Check if wave is missing one arrangement present in crab
  aux_arranget_missing = aux_crab_freqs_df[!aux_crab_freqs_df$arrangement %in% aux_wave_freqs_df$arrangement, 'arrangement']
  if(length(aux_arranget_missing) > 0){
    aux_wave_freqs_df = rbind(aux_crab_freqs_df,
                              data.frame(Population = 'Wave',
                                         arrangement = aux_arranget_missing,
                                         count = 0,
                                         frequency = 0))
  }
  
  wave_freqs_df = data.frame(
    Year = aux_wave_freqs_df$Population,
    Frequency = aux_wave_freqs_df$frequency,
    Allele = c(aux_name_top_1, aux_name_top_2),
    Population = 'Wave',
    Color = aux_wave_freqs_df$arrangement
  )
  
  crab_freqs_df = data.frame(
    Year = aux_crab_freqs_df$Population,
    Frequency = aux_crab_freqs_df$frequency,
    Allele = c(aux_name_top_1, aux_name_top_2),
    Population = 'Wave',
    Color = aux_crab_freqs_df$arrangement
  )
  
  df_arrangements_frequencies_wavecrab = rbind(wave_freqs_df, crab_freqs_df)
  
  df_arrangements_frequencies_skerry$NumSnps = num_snps_in_inversion
  df_arrangements_frequencies_wavecrab$NumSnps = num_snps_in_inversion
  
  list_frequencies = list('Skerry' = df_arrangements_frequencies_skerry, 'WaveCrab' = df_arrangements_frequencies_wavecrab)
  
  return(list_frequencies)
}



#########################################################################
# CREATES A DATAFRAME WITH THE FREQUENCY OF THE ARRANGEMENT
# OF THE MOST COMMON HOMOKARYOTYPE IN THE WAVE ECOTYPE, IN THE DIFFERENT
# YEARS OF THE SKERRY INCLUDING CRAB 1992
#########################################################################
generate_simple_inversion_arrangement_tajectory <- function(pca_coord_2pc, num_snps_in_inversion){
  
  #Create the dataframe with the frequency of the most common arangement in crab and wave
  df_arrangements_frequencies <- data.frame(matrix(ncol = 3, nrow = 0))
  
  frequency_df <- pca_coord_2pc %>% 
    group_by(Population, NewGroup) %>% 
    summarise(count = n())
  
  # Identify the kariotype with the largest frequency in wave
  #tmp_wave_data <- frequency_df[frequency_df$Population == 'Wave',]
  tmp_wave_data <- frequency_df[frequency_df$Population == 'Wave',]
  tmp_crab_data <- frequency_df[frequency_df$Population == 'Crab',]
  tmp_wave_arrf_data <- estimate_simple_inversion_arrangement_frequencies(tmp_wave_data)
  tmp_crab_arrf_data <- estimate_simple_inversion_arrangement_frequencies(tmp_crab_data)
  
  tmp_max_wave_data <- tmp_wave_arrf_data[which.max(tmp_wave_arrf_data$frequency),]
  
  #Find the arrangement that is more common in wave than in crab
  #E.g. if the most common arrangement in wave is 1 but it is also the
  #most common arrangement in crab, 
  #then switch the most common arrangement in wave from 1 to 3
  tmp_max_wave_group = tmp_max_wave_data$arrangement
  tmp_crab_max_wave = tmp_crab_arrf_data[tmp_crab_arrf_data$arrangement == tmp_max_wave_group, ]
  if(tmp_crab_max_wave$frequency > tmp_max_wave_data$frequency){
    if(tmp_max_wave_group == 1){
      tmp_max_wave_group = 3
    }else{
      tmp_max_wave_group= 1
    }
  }
  
  #Estimate the frequency of the crab and wave inversion in skerry for different years
  tmp_crab_1992_data <- frequency_df[frequency_df$Population == 'Crab 1992',]
  tmp_skerry_2005_data <- frequency_df[frequency_df$Population == 'Skerry 2005',]
  tmp_skerry_2018_data <- frequency_df[frequency_df$Population == 'Skerry 2018',]
  tmp_skerry_2021_data <- frequency_df[frequency_df$Population == 'Skerry 2021',]
  
  
  #Estimate the frequency of the two arrangements of a simple inversion (1 or 0)
  tmp_crab_1992_arrf_data   <- estimate_simple_inversion_arrangement_frequencies(tmp_crab_1992_data)
  tmp_skerry_2005_arrf_data <- estimate_simple_inversion_arrangement_frequencies(tmp_skerry_2005_data)
  tmp_skerry_2018_arrf_data <- estimate_simple_inversion_arrangement_frequencies(tmp_skerry_2018_data)
  tmp_skerry_2021_arrf_data <- estimate_simple_inversion_arrangement_frequencies(tmp_skerry_2021_data)
  
  
  #Extract the frequency info of the most common allele (arrangement) in the wave ecotype
  tmp_f_1992_wave_alle = tmp_crab_1992_arrf_data  [tmp_crab_1992_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  tmp_f_2005_wave_alle = tmp_skerry_2005_arrf_data[tmp_skerry_2005_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  tmp_f_2018_wave_alle = tmp_skerry_2018_arrf_data[tmp_skerry_2018_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  tmp_f_2021_wave_alle = tmp_skerry_2021_arrf_data[tmp_skerry_2021_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  tmp_f_wave_wave_alle = tmp_wave_arrf_data[tmp_wave_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  tmp_f_crab_wave_alle = tmp_crab_arrf_data[tmp_crab_arrf_data$arrangement  == tmp_max_wave_group,]$frequency
  
  tmp_row_1992_wave_allele = c('C 1992',tmp_f_1992_wave_alle, 'Wave', tmp_max_wave_group, num_snps_in_inversion)
  tmp_row_2005_wave_allele = c('S 2005',tmp_f_2005_wave_alle, 'Wave', tmp_max_wave_group, num_snps_in_inversion)
  tmp_row_2018_wave_allele = c('S 2018',tmp_f_2018_wave_alle, 'Wave', tmp_max_wave_group, num_snps_in_inversion)
  tmp_row_2021_wave_allele = c('S 2021',tmp_f_2021_wave_alle, 'Wave', tmp_max_wave_group, num_snps_in_inversion)
  tmp_row_wave_wave_allele = c('Wave',tmp_f_wave_wave_alle, 'Wave',   tmp_max_wave_group, num_snps_in_inversion)
  tmp_row_crab_wave_allele = c('Crab',tmp_f_crab_wave_alle, 'Wave',   tmp_max_wave_group, num_snps_in_inversion)
  
  df_arrangements_frequencies_skerry = rbind(df_arrangements_frequencies,
                                             tmp_row_1992_wave_allele,
                                             tmp_row_2005_wave_allele,
                                             tmp_row_2018_wave_allele,
                                             tmp_row_2021_wave_allele)
  
  df_arrangements_frequencies_wavecrab = rbind(df_arrangements_frequencies,
                                               tmp_row_wave_wave_allele,
                                               tmp_row_crab_wave_allele) 
  
  
  colnames(df_arrangements_frequencies_skerry) <- c('Year','Frequency', 'Population', 'Color', 'NumSnps')
  colnames(df_arrangements_frequencies_wavecrab) <- c('Year','Frequency', 'Population', 'Color', 'NumSnps')
  
  list_frequencies = list('Skerry' = df_arrangements_frequencies_skerry, 'WaveCrab' = df_arrangements_frequencies_wavecrab)
  
  return (list_frequencies)
}


############################################
# Generates a PCA based on a SNP dataset 
###########################################
generate_pca_inversion <- function(dataset){
  
  #Remove both samples and SNPs with more than 20% missing data
  dataset = remove_incomplete_samples(dataset, 0.2)
  dataset = remove_incomplete_snps(dataset, 0.2)
  
  #Translatete the genotypes in the dataset from e.g. 0/0 to 1,2 or 3.
  dataset_translated = translate_genotypes(dataset)
  
  if(isTRUE(is_replace_na_genotypes)){
    print('Replace NA Genotypes: YES, then impute for NAs the most common genotype')
    dataset_translated = replace_na_genotypes(dataset, dataset_translated)
  }
  
  #Transpose the dataset matrix, so the individuals will be the rows and SNPs (variables) in the columns
  dataset_translated_transposed = t(dataset_translated[,2:dim(dataset_translated)[2]])
  dataset_translated_transposed = matrix(as.numeric(dataset_translated_transposed), ncol = ncol(dataset_translated_transposed))
  snp_names = dataset_translated[, 1]
  colnames(dataset_translated_transposed) <- snp_names
  
  #Perform a PCA as in the FactoMineR package
  res.pca <- PCA(dataset_translated_transposed, graph = FALSE)
  
  #Add the the names of the samples to the results of the PCA
  res.pca[['ind_names']] <- colnames(dataset)[2:dim(dataset)[2]]
  
  return(res.pca)
}


#########################################################
#Remove columns (Snails) with  > 20% of NAs
#########################################################
remove_incomplete_samples <- function(dataset, p_max_na_snps){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  cols_keep = c(1)
  n_dataset = c()
  #Samples must have info for > 80% of the SNPs
  threshold_na = trunc(nrows * p_max_na_snps) #Estimate the threshold of allowed NAs
  for (c in c(2:ncols)){
    col = dataset[,c]
    na_cols = col[col == "NA"]
    #Remove snails with >= 20% of SNP with NA genotypes
    if(length(na_cols) <= threshold_na){
      cols_keep = c(cols_keep, c)
    }else{
      print(paste("Col ", colnames(dataset)[c], "eliminated, # of NA: ", length(na_cols), "/", nrows))
    }
  }
  n_dataset = dataset[, cols_keep]
  return(n_dataset)
}

#######################################################################
#Remove rows (SNPs) with  > 20% of snails with NA genotypes
#######################################################################
remove_incomplete_snps <- function(dataset, p_max_na_snails){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  n_samples = ncols - 1
  threshold_na = trunc(n_samples * p_max_na_snails) #Estimate the threshold of allowed NA snails per SNP
  rows_keep = c()
  n_dataset = c()
  
  #Create a list of thresholds of missing snails within each population
  populations = unlist(strsplit(sample_names_analysis,"\\|"))
  threshold_na_populations = c()
  col_numbers_populations = list()
  tmp_colnames_analysis = colnames(dataset)
  for(p in c(1:length(populations))){
    tmp_l_cols_pop = grep(populations[p], tmp_colnames_analysis)
    col_numbers_populations[[p]] = tmp_l_cols_pop
    #Allow as much as e.g. 20% of NA snails per SNP within one population
    tmp_threshold_pop = trunc(length(tmp_l_cols_pop) * p_max_na_snails)
    threshold_na_populations = c(threshold_na_populations, tmp_threshold_pop)
  }
  
  for (i in c(1:nrows)){
    row = dataset[i,]
    tmp_na_row = row[row == 'NA']
    keep_SNP = TRUE
    #Examine if there are more than p_max_na_snails% of NA snails for the SNP in total
    if(length(tmp_na_row) > threshold_na){
      keep_SNP = FALSE
      print(n_samples)
      print(paste("->SNP ", i, "eliminated, # of NA Snails: ", length(tmp_na_row),"/",n_samples))
    }else{
      #Examine if there are more than p_max_na_snails% of NA snails 
      # within at least one population of the current time point
      for(p in c(1:length(populations))){
        tmp_l_cols_pop = unlist(col_numbers_populations[[p]])
        tmp_pop_row = row[tmp_l_cols_pop]
        tmp_na_pop = tmp_pop_row[tmp_pop_row == 'NA']
        if(length(tmp_na_pop) > threshold_na_populations[p]){
          keep_SNP = FALSE
          print(paste("----> SNP ", i, "eliminated, # of NA Snails: ", length(tmp_na_row), "/", length(tmp_l_cols_pop), "in POP ", populations[p]))
        }
      }
    }
    
    if(keep_SNP){
      rows_keep = c(rows_keep, i)
    }
  }
  n_dataset = dataset[rows_keep, ]
  return(n_dataset)
}


######################################################
# Translates 0/0 into 0, 0/1 into 1, and 1/1 into 2
######################################################
translate_genotypes <- function(dataset){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  n_samples = ncols - 1
  n_dataset = c()
  
  for (i in c(1:nrows)){
    changed = FALSE
    n_0 = 0
    n_1 = 0
    n_2 = 0
    n_a = 0
    row = dataset[i,]
    #n_row = c(toString(row[1]), rep(0, ncols))
    n_row = rep(0, n_samples)
    for (j in c(2:ncols)){
      genotype = toString(row[j])
      if(genotype == '0/0'){
        n_row[j - 1] = 0
        n_0 = n_0 + 1
      }else{
        if(genotype == '0/1'){
          n_row[j - 1] = 1
          n_1 = n_1 + 1
        }else{
          if(genotype == '1/1'){
            n_row[j - 1] = 2
            n_2 = n_2 + 1
          }else{
            n_row[j - 1] = 3
            n_a = n_a + 1
          }
        }
      }
    }
    
    n_row = append(toString(row[1]), n_row)
    n_dataset = rbind(n_dataset, n_row)
    
  }
  return(n_dataset)
}

######################################################
# Replace/impute NA genotypes (number 3 in the translated) 
# with the most common genotype within the population-year
######################################################
replace_na_genotypes <- function(dataset,dataset_translated){
  
  nrows = dim(dataset_translated)[1]
  tmp_colnames_analysis = colnames(dataset)
  populations = unlist(strsplit(sample_names_analysis,"\\|"))
  col_numbers_populations = list()
  
  # Create a list with the coordinates of columns per population
  for(p in c(1:length(populations))){
    tmp_l_cols_pop = grep(populations[p], tmp_colnames_analysis)
    col_numbers_populations[[p]] = tmp_l_cols_pop
  }
  
  #Examine row by row the NA genotypes within each population
  for (i in c(1:nrows)){
    row = dataset_translated[i,]
    for(p in c(1:length(populations))){
      tmp_l_cols_pop = unlist(col_numbers_populations[[p]])
      tmp_pop_row = row[tmp_l_cols_pop]
      tmp_cols_na_genotypes = grep(3, tmp_pop_row)
      #If there is at least one NA genotype (3), replace
      #the 3 for the most common genotype, either 0, 1 or 2.
      if(length(tmp_cols_na_genotypes) > 0){
        
        tmp_most_common_genotoype = names(sort(table(tmp_pop_row[tmp_pop_row != 3]),decreasing=TRUE))[1]
        #Obtain the coordinates of the genotype 3 in the full dataset
        tmp_nacols_full_dataset = col_numbers_populations[[p]][tmp_cols_na_genotypes]
        dataset_translated[i,tmp_nacols_full_dataset] = tmp_most_common_genotoype
        print(paste("REPLACE:", i, "pop", populations[p], "#NAs:", length(tmp_cols_na_genotypes), "Replaced by", tmp_most_common_genotoype, "Coord:", paste(tmp_nacols_full_dataset, collapse = "|")))
      }
    }
  }
  
  return (dataset_translated)
}


#############################################
# GENERATES A DATAFRAME WITH THE INDIVIDUALS 
# CLASSIFIED IN GROUPS BASED ON K-MEANS ALGORITHM
# num_components is how many principal components
# must be used to estimate the clusters, for example
# only the first axis
#############################################
generate_kmeans_clustering <- function(res_pca, k, num_components, group_populations){
  
  pca_coord_2pc <- as.data.frame(res_pca$ind$coord[,1:2])
  #set.seed(42)
  
  #In most cases the clustering must be estimated only based on
  #The first PC, for example for simple inversions.
  #For complex inversions, the clusters must be estimated
  #based on the first and second PC. 
  #Therefore we use the parameter num_components
  kmeans_92 = kmeans(pca_coord_2pc[,1:num_components], centers = k, nstart = 20)
  groups <- unname(kmeans_92$cluster)
  centroids <- kmeans_92$centers
  
  rownames(pca_coord_2pc) <- res_pca$ind_names
  if(group_populations){
    pca_pop_names <- translate_group_pop_names(res_pca$ind_names)
  }else{
    pca_pop_names <- translate_pop_names(res_pca$ind_names, FALSE)
  }
  
  pca_coord_2pc <- cbind(pca_pop_names, groups, pca_coord_2pc)
  
  rownames(pca_coord_2pc) <- res_pca$ind_names
  
  colnames(pca_coord_2pc) <- c('Population', 'Group','PC1','PC2')
  
  
  pca_coord_2pc <- pca_coord_2pc[order(pca_coord_2pc$PC1, pca_coord_2pc$Group),]
  
  # Redefine the ID of the groups (clusters)
  # So, the most left cluster will the cluster 1 instead of the random number
  # assigned by k-means algotithm and the next most left will be group 2 and so on.
  # this is to keep a standard in coloring the groups
  tmp_unique_groups = unique(pca_coord_2pc$Group)
  aux_col_new_group = rep(0,length(pca_coord_2pc$Group))
  aux_col_centroid_x = rep(0,length(pca_coord_2pc$Group))
  aux_col_centroid_y = rep(0,length(pca_coord_2pc$Group))
  
  for(i in seq(1,length(tmp_unique_groups))){
    tmp_replace_rows = grep(tmp_unique_groups[i],pca_coord_2pc$Group)
    aux_col_new_group[tmp_replace_rows] = i
    
    #We need to keep the info of the centroids 
    #to find the most extreme clusters (homokaryotypes) in complex inversion
    tmp_row_centroids = kmeans_92$centers[tmp_unique_groups[i],]
    aux_col_centroid_x[tmp_replace_rows] = tmp_row_centroids[1]
    aux_col_centroid_y[tmp_replace_rows] = tmp_row_centroids[2]
  }
  
  pca_coord_2pc <- cbind(pca_coord_2pc, NewGroup = aux_col_new_group, Centroid_x = aux_col_centroid_x, Centroid_y = aux_col_centroid_y)
  
  return(pca_coord_2pc)
}


######################################################
# Translate the name of the populations from the name of the SNPs
# separate_reference_pops, when TRUE, add the year to the Wave and Crab 2018 and 2021
# when FALSE, Wave and Crab (2018 & 2021) will be reported as just Wave or just Crab
######################################################
translate_pop_names <- function(sample_names, separate_reference_pops){
  pop_names = c()
  for(s_name in sample_names){
    tmp_name = ""
    if(grepl("SKR-05|SRK-21|ExpSkerry", s_name)){
      tmp_name = "Skerry"
      if(grepl("SKR-05", s_name))
        tmp_name = paste(tmp_name, "2005")
      if(grepl("SRK-21", s_name))
        tmp_name = paste(tmp_name, "2021")
      if(grepl("ExpSkerry", s_name))
        tmp_name = paste(tmp_name, "2018")
      
    }else{
      if(grepl("DO-21|DO-92|DonCrab", s_name)){
        tmp_name = "Crab"
        if(grepl("DO-92", s_name))
          tmp_name = paste(tmp_name, "1992")
        if(separate_reference_pops){
          if(grepl("DO-21", s_name))
            tmp_name = paste(tmp_name, "2021")
          if(grepl("DonCrab", s_name))
            tmp_name = paste(tmp_name, "2018")
        }
      }else{
        if(grepl("W-REF|RefW", s_name)){
          tmp_name =  "Wave"
          if(separate_reference_pops){
            if(grepl("W-REF", s_name))
              tmp_name = paste(tmp_name, "2021")
            if(grepl("RefW", s_name))
              tmp_name = paste(tmp_name, "2018")
          }
        }else{
          print("OTRO")
          print(s_name)
          tmp_name =  "Others"
        }
      }
    }
    pop_names = c(pop_names, tmp_name)
  }
  return(pop_names)
}


######################################################
# Translate the name of the group of the populations 
# from the name of the SNPs
######################################################
translate_group_pop_names <- function(sample_names){
  pop_names = c()
  for(s_name in sample_names){
    tmp_name = ""
    if(grepl("SKR-05|SRK-21|ExpSkerry", s_name)){
      tmp_name = "Skerry"
    }else{
      if(grepl("DO-21|DO-92|DonCrab", s_name)){
        tmp_name = "Crab"
      }else{
        if(grepl("W-REF|RefW", s_name)){
          tmp_name =  "Wave"
        }else{
          print("OTRO")
          print(s_name)
          tmp_name =  "Others"
        }
      }
    }
    pop_names = c(pop_names, tmp_name)
  }
  return(pop_names)
}


#The heterokaryotype is the centroid that has the
#shorter distance to the reference homokaryotype. 
#So we need to estimate the length of the segments between
#two candidates and the reference homokaryotype.
#This function returns the id of the group of the heterokaryotype
determine_heterokaryotype <- function(p_haplotype_data, p_heterok_pool, p_ngroup_homok, p_group_candidate_1, p_group_candidate_2){
  
  coord_ref_homok <- unique(p_haplotype_data[p_haplotype_data$NewGroup == p_ngroup_homok,c("Centroid_x","Centroid_y")])
  coord_het1 <- unique(p_haplotype_data[p_haplotype_data$NewGroup == p_group_candidate_1,c("Centroid_x","Centroid_y")])
  coord_het2 <- unique(p_haplotype_data[p_haplotype_data$NewGroup == p_group_candidate_2,c("Centroid_x","Centroid_y")])
  
  length_C_het1 <- sqrt(sum((coord_ref_homok - coord_het1)^2))
  length_C_het2 <- sqrt(sum((coord_ref_homok - coord_het2)^2))
  
  #What is the shortest segment? that is the heterokariotype CB
  #Therefore, the heterokariotype AB is the other group
  #print(length_C_het1)
  #print(length_C_het2)
  aux_group_het = p_group_candidate_1
  if(length_C_het2 < length_C_het1){
    aux_group_het = p_group_candidate_2
  }
  
  return(aux_group_het)
}


##############################################################################################################
# Based on the coordinates of the centroids of the k-means clusters,
# identify three homokariotypes (extremes in the triangle)
# and three heterokariotypes (edges of the triangle).
# This function returns a list key value with the number of the cluster of the kariotypes
##############################################################################################################
determine_complex_karyotypes <- function(p_haplotype_data){
  #Find the minimum and maximum values of the x coordinates
  #that will help to identify two homokariotypes (A and B)
  min_centroid_xval <- min(p_haplotype_data$Centroid_x)
  max_centroid_xval <- max(p_haplotype_data$Centroid_x)
  
  ngroup_homok_A <- max(p_haplotype_data[p_haplotype_data$Centroid_x == min_centroid_xval,"NewGroup"])
  ngroup_homok_B <- max(p_haplotype_data[p_haplotype_data$Centroid_x == max_centroid_xval,"NewGroup"])
  
  
  #Find the minimum and maximum values of the y coordinates
  #that are not the homokaryotype A or B.
  #It will help to identify the third homokaryotype 
  min_centroid_yval <- min(p_haplotype_data[p_haplotype_data$NewGroup != ngroup_homok_A 
                                            & p_haplotype_data$NewGroup != ngroup_homok_B,
                                            "Centroid_y"])
  max_centroid_yval <- max(p_haplotype_data[p_haplotype_data$NewGroup != ngroup_homok_A 
                                            & p_haplotype_data$NewGroup != ngroup_homok_B,
                                            "Centroid_y"])
  
  #The centroid of the homokaryotype C is the largest value along the Y axis
  centroid_yval <- max_centroid_yval
  if(abs(min_centroid_yval) > abs(max_centroid_yval)){
    centroid_yval <- min_centroid_yval
  }
  
  #Find the name of the group of the homokaryotype C that is the group
  #with the largest or shortest centroid along Y axis, but also that is not Homok A nor B.
  ngroup_homok_C <- max(p_haplotype_data[p_haplotype_data$Centroid_y == centroid_yval 
                                         & p_haplotype_data$NewGroup != ngroup_homok_A 
                                         & p_haplotype_data$NewGroup != ngroup_homok_B,
                                         "NewGroup"])
  
  
  
  
  #Identify the group of the heterokaryotypes
  heterok_pool <- unique(p_haplotype_data[p_haplotype_data$NewGroup != ngroup_homok_A 
                                          & p_haplotype_data$NewGroup != ngroup_homok_B
                                          & p_haplotype_data$NewGroup != ngroup_homok_C,
                                          c("NewGroup", "Centroid_x", "Centroid_y")])
  
  heterok_pool <- heterok_pool[order(heterok_pool$Centroid_x),]
  
  ngroup_hetAC = determine_heterokaryotype(p_haplotype_data, heterok_pool, ngroup_homok_C, heterok_pool$NewGroup[1], heterok_pool$NewGroup[2])
  ngroup_hetBC = determine_heterokaryotype(p_haplotype_data, heterok_pool, ngroup_homok_C, heterok_pool$NewGroup[2], heterok_pool$NewGroup[3])
  #aux_group_hetAB = determine_heterokaryotype(p_haplotype_data, heterok_pool, ngroup_homok_B, heterok_pool$NewGroup[1], heterok_pool$NewGroup[2])
  ngroup_hetAB = heterok_pool[heterok_pool$NewGroup != ngroup_hetAC & heterok_pool$NewGroup != ngroup_hetBC, "NewGroup"]
  
  # create the list of kariotype groups
  list_karyotypes <- list()
  
  # Build up key value pairs of Homokaryotypes
  list_karyotypes[[ "A" ]] <- ngroup_homok_A
  list_karyotypes[[ "B" ]] <- ngroup_homok_B
  list_karyotypes[[ "C" ]] <- ngroup_homok_C
  
  # Build up key value pairs of Heterokaryotypes
  list_karyotypes[[ "AC" ]] <- ngroup_hetAC
  list_karyotypes[[ "BC" ]] <- ngroup_hetBC
  list_karyotypes[[ "AB" ]] <- ngroup_hetAB
  
  return (list_karyotypes)
}
