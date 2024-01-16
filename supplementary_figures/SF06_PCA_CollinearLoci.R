########################## PCA Collinear Loci ############################
#Performns a PCA on both control and spatial outlier SNPs and 
#generates a PDF file with the scatter plots of the clustering pattern
#
#author: "Diego Garcia"
#date: "2023-02-09"
################################################################################
  
#Set the working directory
setwd('.')

#Load libraries
library(ggplot2)
library(ggfortify)
library("FactoMineR")
library("factoextra")
library("cowplot")
library(ggridges)
library(RColorBrewer)
library('stringr')
library(dplyr)

output_filename = 'PCA_Swedish_loci_Final' #Name of the output PDF file with the plot

#Define the path and name of the files that contain control and spatial outliers SNP genotypes
swedish_neutral_snp_file <- "../data/SkerryExperiment_Neutral_NOLG12_Swedenfilter.txt"
swedish_outlier_snp_file <- "../data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt"

#Read and load the data files with the genotypes
neutral <- read.table(swedish_neutral_snp_file, header = T, check.names = FALSE)
outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE)

# Define the colours that the PCA scatter plot will use for each population
colors_crab_skerry_wave = c("Skerry 2005" =		  "#DE4DC1",
                            "Skerry 2018" =		  "#83D873",
                            "Skerry 2021" =		  "#9266C1",
                            "Donor Crab 1992" =	"#f60000",
                            "Donor Crab 2018" =	"#e69138",
                            "Donor Crab 2021" =	"#c90076",
                            "Reff. Wave 2018" =	"#0068ff",
                            "Reff. Wave 2021" =	"#00c7ff",
                            "Donor Crab" =		  "#E09157",
                            "Reff. Wave" =		  "#73BADA")

#neutral = original_neutral[original_neutral$Contig %in% swedish_neutral$cat, ]
#outliers = original_outliers[original_outliers$Contig %in% swedish_outliers$cat, ]

l_cols = colnames(neutral)

#Define the keywords of the sample names (columns) that will be used in the analysis
sample_names_1992 = "DO.92"
sample_names_2005 = "SKR.05"
sample_names_2018 = "DonCrab|ExpSkerry|RefW"
sample_names_2021 = "DO.21|SRK.21|W.REF"

#Merge the column names that will be used in the analysis, this pattern
#is used to subset the columns
sample_names_analysis = paste(sample_names_1992, sample_names_2005, sample_names_2018,sample_names_2021,sep = "|")
# sample_names_analysis = paste(sample_names_1992, sample_names_2005, sample_names_2021,sep = "|")
# sample_names_analysis = "DO-92|DonCrab|RefW|DO-21|SKR-05|W-REF"
l_cols_pca = grep(sample_names_analysis, l_cols)
l_cols_pca = c( grep('cat',l_cols), l_cols_pca)

pca_samples_neutral = neutral[,l_cols_pca]
pca_samples_outlier = outliers[,l_cols_pca]
# PCA ONLY ON THE CANDIDATE LOCI
#pca_samples_outlier = outliers[match(tmp_colnames_outliers, outliers$Contig), ]

max_na_snails = 0.05 #SNPs must have info for > 95% of the snails
max_na_snps = 0.05 #Snails must have info for > 95% of the SNPs
replace_na_genotypes_flag = TRUE #Imputes the most common genotype for a missing SNP data

#Clean neutral and outliers by removing samples with more than 5% of NA genotypes
pca_samples_neutral_clean = remove_incomplete_samples(pca_samples_neutral, max_na_snps)
pca_samples_neutral_clean = remove_incomplete_snps(pca_samples_neutral_clean, max_na_snails)

pca_samples_outlier_clean = remove_incomplete_samples(pca_samples_outlier, max_na_snps)
pca_samples_outlier_clean = remove_incomplete_snps(pca_samples_outlier_clean, max_na_snails)

#Invoke the function that performs a PCA on the given dataset
neutral.res.pca <- generate_pca(pca_samples_neutral_clean)
outliers.res.pca <- generate_pca(pca_samples_outlier_clean)


split_ref_populations = FALSE #When FALSE, Crab and Wave populations will be combined

#Generate the scatter plots of the PCAs
pltn <- plot_pca_figures(neutral.res.pca, pca_samples_neutral_clean, "Control", split_ref_populations)
plto <- plot_pca_figures(outliers.res.pca, pca_samples_outlier_clean, "Spatial Outlies", split_ref_populations)

######## GENERATE PDF OUTPUT #########################################
pdf(paste0(output_filename,'.pdf'),width = 10, height = 5)

#Create grid with two cells: neutral and outliers PCA
plt_row <- plot_grid(
  pltn + theme(legend.position="none"),
  plto + theme(legend.position="none"),
  nrow = 1,
  ncol = 2
)

# extract a legend that is laid out horizontally
legend_b <- get_legend(
  pltn + 
    theme(legend.position = "bottom")
)

# Plot the figures on the top row
# and the legend on the bottom row
plot_grid(
  plt_row,
  legend_b,
  ncol = 1,
  rel_heights = c(1, .2))
dev.off()


##########################################################
#Remove columns (Snails) with  >= 5% of NAs
##########################################################
remove_incomplete_samples <- function(dataset, p_max_na_snps){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  cols_keep = c(1)
  n_dataset = c()
  #Samples must have info for > 95% of the SNPs
  threshold_na = trunc(nrows * p_max_na_snps) #Allow as much as 10% of NA SNPs per sample
  for (c in c(2:ncols)){
    col = dataset[,c]
    na_cols = col[col == "NA"]
    #Remove snails with > 5% of SNP with NA genotypes
    if(length(na_cols) <= threshold_na){
      cols_keep = c(cols_keep, c)
    }else{
      print(paste("Col ", colnames(dataset)[c], "eliminated, # of NA: ", length(na_cols), "/", nrows))
    }
  }
  n_dataset = dataset[, cols_keep]
  return(n_dataset)
}

##################################################################
#Remove rows (SNPs) with  >= 5% of snails with NA genotypes
##################################################################
remove_incomplete_snps <- function(dataset, p_max_na_snails){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  n_samples = ncols - 1
  threshold_na = trunc(n_samples * p_max_na_snails) #Allow as much as 5% of NA snails per SNP
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
    #Allow as much as 5% of NA snails per SNP within one population
    tmp_threshold_pop = trunc(length(tmp_l_cols_pop) * p_max_na_snails)
    threshold_na_populations = c(threshold_na_populations, tmp_threshold_pop)
  }
  
  
  for (i in c(1:nrows)){
    row = dataset[i,]
    
    tmp_na_row = row[row == 'NA']
    keep_SNP = TRUE
    #Examine if there are more than 5% of NA snails for the SNP in total
    if(length(tmp_na_row) > threshold_na){
      keep_SNP = FALSE
      print(n_samples)
      print(paste("->SNP ", i, "eliminated, # of NA Snails: ", length(tmp_na_row),"/",n_samples))
    }else{
      #Examine if there are more than 5% of NA snails 
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
    
    # PRPCOMP REQUIERES TO ELIMINATE ROWS THAT DO NOT HAVE VARIATION
    # if(keep_SNP){
    #   #Exclude SNPs with no variation accros all samples since variance is 0 and it is
    #   #a scaling problem before ploting the PCA
    #   num_unique_elements = length(unique(unlist(row, use.names = FALSE)))
    #   if(num_unique_elements > 2){
    #     rows_keep = c(rows_keep, i)
    #   }else{
    #     print(paste("----> SNP ", i, "eliminated, there is no variation accross the samples"))
    #   }
    # }
    
    
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
# Teplace NA genotypes (number 3 in the translated) 
# with the most common genotype within the population
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



######################################################
# Translate the name of the populations from the name of the SNPs
######################################################
translate_pop_names <- function(sample_names, p_split_source_populations){
  pop_names = c()
  for(s_name in sample_names){
    tmp_name = ""
    if(grepl("SKR-05|SRK-21|ExpSkerry|SKR.05|SRK.21", s_name)){
      tmp_name = "Skerry"
      if(grepl("SKR-05|SKR.05", s_name))
        tmp_name = paste(tmp_name, "2005")
      if(grepl("SRK-21|SRK.21", s_name))
        tmp_name = paste(tmp_name, "2021")
      if(grepl("ExpSkerry", s_name))
        tmp_name = paste(tmp_name, "2018")
      
    }else{
      if(grepl("DO-21|DO-92|DonCrab|DO.21|DO.92", s_name)){
        tmp_name = "Donor Crab"
        if(p_split_source_populations){
          if(grepl("DO-21|DO.21", s_name))
            tmp_name = paste(tmp_name, "2021")
          if(grepl("DO-92|DO.92", s_name))
            tmp_name = paste(tmp_name, "1992")
          if(grepl("DonCrab", s_name))
            tmp_name = paste(tmp_name, "2018")
        }
      }else{
        if(grepl("W-REF|RefW|W.REF", s_name)){
          tmp_name =  "Reff. Wave"
          if(p_split_source_populations){
            if(grepl("W-REF|W.REF", s_name))
              tmp_name = paste(tmp_name, "2021")
            if(grepl("RefW", s_name))
              tmp_name = paste(tmp_name, "2018")
          }
        }else{
          print("OTHER")
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
# Generate the PCA on a dataset
######################################################
generate_pca <- function(dataset){
  dataset_translated = translate_genotypes(dataset)
  if(isTRUE(replace_na_genotypes_flag)){
    dataset_translated = replace_na_genotypes(dataset, dataset_translated)
  }
  dataset_translated_transposed = t(dataset_translated[,2:dim(dataset_translated)[2]])
  dataset_translated_transposed = matrix(as.numeric(dataset_translated_transposed), ncol = ncol(dataset_translated_transposed))
  snp_names = dataset_translated[, 1]
  colnames(dataset_translated_transposed) <- snp_names
  
  res.pca <- PCA(dataset_translated_transposed, graph = FALSE, ncp = 200)
  return(res.pca)
}



######################################################
# Generate a scatter plot with icons for each group of snails
######################################################
plot_pca_figures <- function(res.pca, dataset, title, split_ref_populations){
  
  pop_names = translate_pop_names(colnames(dataset)[2:dim(dataset)[2]], split_ref_populations)
  variance_pc1 <- res.pca$eig[1,2]
  variance_pc2 <- res.pca$eig[2,2]
  
  plt <- fviz_pca_ind(res.pca,
                      geom.ind = "point", # show points only (nbut not "text")
                      col.ind = factor(pop_names),
                      fill.ind = factor(pop_names),
                      pointsize = 1,
                      alpha.ind = 0.8,
                      mean.point = FALSE,
                      addEllipses = TRUE, # Concentration ellipses
                      ellipse.alpha = 0.1,
                      #ellipse.level = 0.975,
                      repel = TRUE,
                      title = title,
                      ggtheme = theme_minimal(),
                      axes.linetype=NA)+
    scale_shape_manual(values=c(1,2,3,4,5,6,7,8,9))+
    scale_color_manual(values = colors_crab_skerry_wave) +
    scale_fill_manual(values = colors_crab_skerry_wave) +
    ylim(-11,11)+
    xlim(-18,18)+
    theme(
      legend.title = element_blank(),
      plot.title = element_text(size=10, face = "bold", hjust = 0.5),
      axis.title.x=element_text(size=8, face = "bold"),
      axis.title.y=element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=0.1))+
    labs(x = paste('PC1 (', round(variance_pc1, digits=1), '%)', sep=""), 
         y = paste('PC2 (', round(variance_pc2, digits=1), '%)', sep=""))
  
  return(plt)
}