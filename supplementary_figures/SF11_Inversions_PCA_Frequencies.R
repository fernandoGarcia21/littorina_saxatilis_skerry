############## PCA ON INVERSIONS + INVERSION FREQUENCIES #######################
#Performns a PCA on inversion SNPs, identify the clusters by k-means,
#estimate the frequencies of the Wave arrangement and 
#Crab arrangement (only for complex inversions) based on the clustering pattern,
#and generates a PDF file with the scatter plots of the clustering pattern
#as well as a barplot with the frequencies
#
#author: "Diego Garcia"
#date: "05-Oct-2022"
################################################################################


#Set the working directory
setwd('.')

#Load libraries
library("ggplot2")
library('stringr')
library("FactoMineR")
library("factoextra")
library("cowplot")
library("dplyr")
library(forcats)

#Define the path and name of the files that contain inversion SNP genotypes and coordinates
inv_file <- "../Data/Inv_Genotypes.txt" #SNP genotypes of the inversions
inv_coordinates <- "../Data/Inv_Coordinates.txt" #Coordinates of the inversions (e.g. start and end in the genetic map)
contig_pos_gm_file <- "../Data/contigPositionsGM.txt" #Coordinates of the contigs in the genetic map

#Read files and load data
inv_data <- read.table(inv_file, header = T, check.names = FALSE)
inv_coordinates_data <- read.table(inv_coordinates, header = T, check.names = FALSE)
pos_gm_data <- read.table(contig_pos_gm_file, header = T, check.names = FALSE)

l_cols = colnames(inv_data)

#Define the colour scheme of the scatter plots and bar plots
palette_all <- c('#85C9E8','#EF9797', '#83D873',
                 '#9676CC', '#729A6B', '#BF56AE')

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
is_replace_na_genotypes = TRUE #Impute the NA genotypes with the most common genotype in the population
pca_allp_lg_1_1 <- generate_pca_inversion(tmp_lg_data_1_1)
pca_allp_lg_1_2 <- generate_pca_inversion(tmp_lg_data_1_2)
pca_allp_lg_2 <- generate_pca_inversion(tmp_lg_data_2)
pca_allp_lg_4 <- generate_pca_inversion(tmp_lg_data_4)
#pca_allp_lg_5 <- generate_pca_inversion(tmp_lg_data_5) #Excluded putative inversion
pca_allp_lg_6 <- generate_pca_inversion(tmp_lg_data_6)
pca_allp_lg_7_1 <- generate_pca_inversion(tmp_lg_data_7_1)
pca_allp_lg_7_2 <- generate_pca_inversion(tmp_lg_data_7_2)
pca_allp_lg_9 <- generate_pca_inversion(tmp_lg_data_9)
pca_allp_lg_10_1 <- generate_pca_inversion(tmp_lg_data_10_1)
pca_allp_lg_10_2 <- generate_pca_inversion(tmp_lg_data_10_2)
pca_allp_lg_11 <- generate_pca_inversion(tmp_lg_data_11)
pca_allp_lg_14_1 <- generate_pca_inversion(tmp_lg_data_14_1)
#pca_allp_lg_14_2 <- generate_pca_inversion(tmp_lg_data_14_2) #Excluded putative inversion
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
merge_skerry_years = TRUE 

coord_2pc_k_1_1  <- generate_kmeans_clustering(pca_allp_lg_1_1,   3, 1, merge_skerry_years)
coord_2pc_k_1_2  <- generate_kmeans_clustering(pca_allp_lg_1_2,   3, 1, merge_skerry_years)
coord_2pc_k_2    <- generate_kmeans_clustering(pca_allp_lg_2,     3, 1, merge_skerry_years)
coord_2pc_k_4    <- generate_kmeans_clustering(pca_allp_lg_4,     3, 1, merge_skerry_years)
#coord_2pc_k_5    <- generate_kmeans_clustering(pca_allp_lg_5,    3, 1, merge_skerry_years) #excluded putative inversion
coord_2pc_k_6    <- generate_kmeans_clustering(pca_allp_lg_6,     6, 2, merge_skerry_years)
coord_2pc_k_7_1  <- generate_kmeans_clustering(pca_allp_lg_7_1,   3, 1, merge_skerry_years)
coord_2pc_k_7_2  <- generate_kmeans_clustering(pca_allp_lg_7_2,   3, 1, merge_skerry_years)
coord_2pc_k_9    <- generate_kmeans_clustering(pca_allp_lg_9,     3, 1, merge_skerry_years)
coord_2pc_k_10_1 <- generate_kmeans_clustering(pca_allp_lg_10_1,  3, 1, merge_skerry_years)
coord_2pc_k_10_2 <- generate_kmeans_clustering(pca_allp_lg_10_2,  3, 1, merge_skerry_years)
coord_2pc_k_11   <- generate_kmeans_clustering(pca_allp_lg_11,    3, 1, merge_skerry_years)
coord_2pc_k_14_1 <- generate_kmeans_clustering(pca_allp_lg_14_1,  6, 2, merge_skerry_years)
#coord_2pc_k_14_2 <- generate_kmeans_clustering(pca_allp_lg_14_2, 3, 1, merge_skerry_years) #excluded putative inversion
coord_2pc_k_17   <- generate_kmeans_clustering(pca_allp_lg_17,    3, 1, merge_skerry_years)

#2. Create the scatter plots of each LG based on the coordinates of the PC1 and PC2
plt_k_1_1  <- plot_pca_by_kmeans_clustering(pca_allp_lg_1_1,   coord_2pc_k_1_1 , 'LGC1.1'   )
plt_k_1_2  <- plot_pca_by_kmeans_clustering(pca_allp_lg_1_2,   coord_2pc_k_1_2 , 'LGC1.2'   )
plt_k_2    <- plot_pca_by_kmeans_clustering(pca_allp_lg_2,     coord_2pc_k_2   , 'LGC2.1'   )
plt_k_4    <- plot_pca_by_kmeans_clustering(pca_allp_lg_4,     coord_2pc_k_4   , 'LGC4.1'   )
plt_k_6    <- plot_pca_by_kmeans_clustering(pca_allp_lg_6,     coord_2pc_k_6   , 'LGC6.1/2' )
plt_k_7_1  <- plot_pca_by_kmeans_clustering(pca_allp_lg_7_1,   coord_2pc_k_7_1 , 'LGC7.1'   )
plt_k_7_2  <- plot_pca_by_kmeans_clustering(pca_allp_lg_7_2,   coord_2pc_k_7_2 , 'LGC7.2'   )
plt_k_9    <- plot_pca_by_kmeans_clustering(pca_allp_lg_9,     coord_2pc_k_9   , 'LGC9.1'   )
plt_k_10_1 <- plot_pca_by_kmeans_clustering(pca_allp_lg_10_1,  coord_2pc_k_10_1, 'LGC10.1'  )
plt_k_10_2 <- plot_pca_by_kmeans_clustering(pca_allp_lg_10_2,  coord_2pc_k_10_2, 'LGC10.2'  )
plt_k_11   <- plot_pca_by_kmeans_clustering(pca_allp_lg_11,    coord_2pc_k_11  , 'LGC11.1'  )
plt_k_14_1 <- plot_pca_by_kmeans_clustering(pca_allp_lg_14_1,  coord_2pc_k_14_1, 'LGC14.1/2')
plt_k_17   <- plot_pca_by_kmeans_clustering(pca_allp_lg_17,    coord_2pc_k_17  , 'LGC17.1'  )


#3. Estimate the frequency of each cluster within each population and generate the bar plot of each LG
plt_k_freq_1_1  <- plot_freq_by_kmeans_clustering(coord_2pc_k_1_1 , FALSE, 'LGC1.1'   )
plt_k_freq_1_2  <- plot_freq_by_kmeans_clustering(coord_2pc_k_1_2 , FALSE, 'LGC1.2'   )
plt_k_freq_2    <- plot_freq_by_kmeans_clustering(coord_2pc_k_2   , FALSE, 'LGC2.1'   )
plt_k_freq_4    <- plot_freq_by_kmeans_clustering(coord_2pc_k_4   , FALSE, 'LGC4.1'   )
plt_k_freq_6    <- plot_freq_by_kmeans_clustering(coord_2pc_k_6   , FALSE, 'LGC6.1/2' )
plt_k_freq_7_1  <- plot_freq_by_kmeans_clustering(coord_2pc_k_7_1 , FALSE, 'LGC7.1'   )
plt_k_freq_7_2  <- plot_freq_by_kmeans_clustering(coord_2pc_k_7_2 , FALSE, 'LGC7.2'   )
plt_k_freq_9    <- plot_freq_by_kmeans_clustering(coord_2pc_k_9   , FALSE, 'LGC9.1'   )
plt_k_freq_10_1 <- plot_freq_by_kmeans_clustering(coord_2pc_k_10_1, FALSE, 'LGC10.1'  )
plt_k_freq_10_2 <- plot_freq_by_kmeans_clustering(coord_2pc_k_10_2, FALSE, 'LGC10.2'  )
plt_k_freq_11   <- plot_freq_by_kmeans_clustering(coord_2pc_k_11  , FALSE, 'LGC11.1'  )
plt_k_freq_14_1 <- plot_freq_by_kmeans_clustering(coord_2pc_k_14_1, FALSE, 'LGC14.1/2')
plt_k_freq_17   <- plot_freq_by_kmeans_clustering(coord_2pc_k_17  , FALSE, 'LGC17.1'  )


################### TO PLOT BAR PLOTS NEXT TO SCATTER PLOTS ###########################
#6.Plot the Scatter plot and the Bar plot one next the other for each LG
pdf('SI_Inversions_PCA_KaryotypeFreqs_Inputed.pdf', width=9, height=6)
#png('Inversions_PCA_KaryotypeFreqs_Final.png', width=9, height=6)
plot_grid(plt_k_1_1,  plt_k_freq_1_1, NULL,  plt_k_1_2,  plt_k_freq_1_2, NULL, plt_k_2,  plt_k_freq_2,
          plt_k_4, plt_k_freq_4, NULL, plt_k_6,  plt_k_freq_6, NULL, plt_k_7_1,  plt_k_freq_7_1,
          plt_k_7_2,  plt_k_freq_7_2, NULL,  plt_k_9,  plt_k_freq_9, NULL, plt_k_10_1, plt_k_freq_10_1,
          plt_k_10_2, plt_k_freq_10_2, NULL, plt_k_11, plt_k_freq_11, NULL, plt_k_14_1, plt_k_freq_14_1, 
          plt_k_17, plt_k_freq_17,
          ncol=8, nrow=5, rel_widths = c(0.19, 0.21, 0.05, 0.19, 0.21, 0.05, 0.19, 0.21))

dev.off()


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


#############################################
# GENERATES A SCATTER PLOT BASED ON A PCA 
# BASED ON K-MEANS ALGORITHM
#############################################
plot_pca_by_kmeans_clustering <- function(res_pca, pca_coord_2pc, title){
  
  variance_pc1 <- res_pca$eig[1,2]
  variance_pc2 <- res_pca$eig[2,2]
  
  plt <- ggplot(pca_coord_2pc, aes(x=PC1, y=PC2))+
    #geom_point(alpha=0.6, size = 1.4, shape=21, color='white', 
    #           stroke = 0.1, aes(fill=factor(Group)))+
    geom_jitter(position=position_jitter(height=0.5, width=0.8),
                alpha=0.9, size = 1.4, shape=21, color='white', 
                stroke = 0.1, aes(fill=factor(NewGroup)))+
    ggtitle(title)+
    labs(x = paste('PC1 (', round(variance_pc1, digits=1), '%)', sep=""), 
         y = paste('PC2 (', round(variance_pc2, digits=1), '%)', sep=""))+
    theme_light() +
    theme(    panel.grid.major = element_blank(),
              legend.position = "none",
              legend.title = element_blank(),
              panel.grid.minor = element_blank(),
              plot.title = element_text(size=6, face = "bold", margin=margin(0,0,2,0)),
              axis.title.x=element_text(size=5, margin=margin(0.5,0,0,0)),
              axis.title.y=element_text(size=5, margin=margin(0,0.5,0,0)),
              axis.text.x=element_text(size=4),
              axis.text.y=element_text(size=4),
              plot.margin = unit(c(0.1,0.1,0.1,0.2), "cm")
    )+
    scale_fill_manual(values = palette_all)
  #scale_colour_manual(values = rev(colors_crab_skerry_wave))
  
  return(plt)
}


#############################################
# GENERATES A BARPLOT PLOT BASED ON THE FREQUENCY 
# OF EACH ARRANGEMENT (KAREOTYPE) WITHIN EACH POPULATION
#############################################
plot_freq_by_kmeans_clustering <- function(pca_coord_2pc, add_title, title){
  
  frequency_df <- pca_coord_2pc %>% 
    group_by(Population, NewGroup) %>% 
    summarise(count = n()) %>% 
    mutate(perc = count/sum(count))
  
  plt <- ggplot(frequency_df, aes(x = factor(Population), y = perc, fill = factor(NewGroup))) +
    geom_bar(stat="identity", width = 0.4, alpha = 0.9, position = position_fill(reverse = TRUE)) +
    coord_flip()+
    labs(y = "Karyotype frequency", x="", fill = "Arrangement") +
    scale_fill_manual(values = palette_all)+
    theme_light()+
    ggtitle("")+
    theme(  panel.grid.major = element_blank(),
            legend.position = "none",
            legend.title = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(size=6, face = "bold", margin=margin(0,0,2,0)),
            axis.title.x=element_text(size=5, margin=margin(0.5,0,0,0)),
            axis.title.y=element_text(size=5, margin=margin(0,0.5,0,0)),
            axis.text.x=element_text(size=4),
            axis.text.y=element_text(size=6),
            plot.margin = unit(c(0.1,0.2,0.1,0.1), "cm"))
  
  if(add_title){
    plt <- plt+ggtitle(title)+theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm"))
  }
  
  return(plt)
  
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