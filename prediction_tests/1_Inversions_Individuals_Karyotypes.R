############## INVERSIONS INDIVIDUALS KARYOTYPES #######################
# Performns a PCA on inversion SNPs, identify the clusters by k-means,
# and determine the karyotype (cluster ID) of each individual. 
# Simple inversions have three karyotypes.
# Complex inversions have six posible karyotypes.
# Write two output files:
# Individuals_karyotypes_inversions.txt is a table with the karyotpes of each individual.
# Complex_Karyotypes_Ids.txt is a translation from 1 to 6 to A,B,C,AC, BC, AB for complex inversions
# These output files are used as input for the Fisher's exact test to validate the hypothesis.
#
# author: "Diego Garcia"
# date: "12-Dec-2023"
################################################################################


#Set the working directory
setwd('.')

#Load libraries
library("FactoMineR")

#Define the path and name of the files that contain inversion SNP genotypes
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


###########################################
# ESTIMATE THE CLUSTERS OF THE PCA BY K-MEANS 
# INSTEAD OF COLORING BY POPULATION 
# AND GENERATE A PCA PLOT AND A BARPLOT FOR EACH INVERSION
###########################################

#1. Perform the clustering on the result of the PCA and return PC1 and PC2 coordinates
# The second parameter is the number of expected clusters, 
# the third parameter is the number of principal components to estimate the clusters.
# Complex inversions need to estimate the clusters from PC1 and PC2

#When True, sampling years are combined by population for Crab, Wave, and Skerry.
#When False, keep the year of the population. e.g. Skerry 2021, Crab 2018
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

#Determine the kariotypes or complex inversions (which are homokaryotypes, which are heterokaryotypes)
kariotypes_k_6 = as.data.frame(determine_complex_karyotypes(coord_2pc_k_6))
kariotypes_k_14_1 = as.data.frame(determine_complex_karyotypes(coord_2pc_k_14_1))
complex_kariotypes = cbind(rbind(kariotypes_k_6,kariotypes_k_14_1),Inversion=c('LGC6.1/2','LGC14.1/2'))

#Print the kariotype id of each individual for each inversion
#these data are used for chi-square/Fisher's, to test the hypothesis accuracy
df_clusters_all_inversions = as.data.frame(matrix(nrow = 0, ncol = 4))
df_clusters_all_inversions = rbind(df_clusters_all_inversions, 
                                   data.frame(Inversion = rep('LGC1.1'   ,nrow(coord_2pc_k_1_1  )), Individual=rownames(coord_2pc_k_1_1  ), Karyotype=coord_2pc_k_1_1  $NewGroup),
                                   data.frame(Inversion = rep('LGC1.2'   ,nrow(coord_2pc_k_1_2  )), Individual=rownames(coord_2pc_k_1_2  ), Karyotype=coord_2pc_k_1_2  $NewGroup),
                                   data.frame(Inversion = rep('LGC2.1'   ,nrow(coord_2pc_k_2    )), Individual=rownames(coord_2pc_k_2    ), Karyotype=coord_2pc_k_2    $NewGroup),
                                   data.frame(Inversion = rep('LGC4.1'   ,nrow(coord_2pc_k_4    )), Individual=rownames(coord_2pc_k_4    ), Karyotype=coord_2pc_k_4    $NewGroup),
                                   data.frame(Inversion = rep('LGC6.1/2' ,nrow(coord_2pc_k_6    )), Individual=rownames(coord_2pc_k_6    ), Karyotype=coord_2pc_k_6    $NewGroup),
                                   data.frame(Inversion = rep('LGC7.1'   ,nrow(coord_2pc_k_7_1  )), Individual=rownames(coord_2pc_k_7_1  ), Karyotype=coord_2pc_k_7_1  $NewGroup),
                                   data.frame(Inversion = rep('LGC7.2'   ,nrow(coord_2pc_k_7_2  )), Individual=rownames(coord_2pc_k_7_2  ), Karyotype=coord_2pc_k_7_2  $NewGroup),
                                   data.frame(Inversion = rep('LGC9.1'   ,nrow(coord_2pc_k_9    )), Individual=rownames(coord_2pc_k_9    ), Karyotype=coord_2pc_k_9    $NewGroup),
                                   data.frame(Inversion = rep('LGC10.1'  ,nrow(coord_2pc_k_10_1 )), Individual=rownames(coord_2pc_k_10_1 ), Karyotype=coord_2pc_k_10_1 $NewGroup),
                                   data.frame(Inversion = rep('LGC10.2'  ,nrow(coord_2pc_k_10_2 )), Individual=rownames(coord_2pc_k_10_2 ), Karyotype=coord_2pc_k_10_2 $NewGroup),
                                   data.frame(Inversion = rep('LGC11.1'  ,nrow(coord_2pc_k_11   )), Individual=rownames(coord_2pc_k_11   ), Karyotype=coord_2pc_k_11   $NewGroup),
                                   data.frame(Inversion = rep('LGC14.1/2',nrow(coord_2pc_k_14_1 )), Individual=rownames(coord_2pc_k_14_1 ), Karyotype=coord_2pc_k_14_1 $NewGroup),
                                   data.frame(Inversion = rep('LGC17.1'  ,nrow(coord_2pc_k_17   )), Individual=rownames(coord_2pc_k_17   ), Karyotype=coord_2pc_k_17   $NewGroup)
                                   
)
#Export the karyotypes ids for complex inversions to a file
write.table(complex_kariotypes, '../Data/Complex_Karyotypes_Ids.txt', quote = FALSE, row.names = FALSE)

#Export the karyotypes of all individuals for all inversions
write.table(df_clusters_all_inversions, '../Data/Individuals_karyotypes_inversions.txt', quote = FALSE, row.names = FALSE)



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
