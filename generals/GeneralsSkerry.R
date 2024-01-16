########################### GENERALS SKERRY################################
# Functions that are invoked from multiple other scripts
#
#author: "Diego Garcia"
#date: "2023-05-11"
############################################################################


######################################################
# Computes the FST for each SNP of two populations 
# and returns an array with FST values
# aaf_pop_A is a vector the alternative allele frequencies of all SNPs in pop A
# aaf_pop_B is a vector the alternative allele frequencies of all SNPs in pop B
# Fst = (Ht-Hs)/Ht
######################################################
compute_fst <- function(aaf_pop_A, aaf_pop_B){
  
  #Compute the frequency of the reference allele for all SNPs
  raf_pop_A = 1 - aaf_pop_A
  raf_pop_B = 1 - aaf_pop_B
  
  #Estimate the expected heterozygosity
  expec_heterozygosity_A = 2 * raf_pop_A * aaf_pop_A
  expec_heterozygosity_B = 2 * raf_pop_B * aaf_pop_B
  
  #Estimate the heterozygosity of the sample (hs)
  h_s = (expec_heterozygosity_A + expec_heterozygosity_B) / 2
  
  #Estimate the total heterozygosity (ht)
  h_t_ref = (raf_pop_A + raf_pop_B) / 2
  h_t_alt = (aaf_pop_A + aaf_pop_B) / 2
  h_t = 2 * h_t_ref * h_t_alt
  #Find the indexes of all ht equal 0
  h_t_0 = which(h_t == 0)
  
  # Fst = (Ht-Hs)/Ht
  fst = (h_t - h_s) / h_t
  print(h_t_0)
  #Replace those Nan values of fst with 0
  fst[h_t_0] = 0
  
  return(fst)
}


######################################################
# Computes the FST for a multiallelic locus
# following the formula of: Fst = (Ht-Hs)/Ht
# where H = 1 - Σ(p^2)
# aaf_pop_A is a vector with the frequencies of two arrangements (out of three) for one complex inversion in pop A
# aaf_pop_B is a vector with the frequencies of two arrangements (out of three) for one complex inversion in pop B
######################################################
compute_fst_multiallelic <- function(aaf_pop_A, aaf_pop_B){
  
  #Add a third element to the arrangement frequencies vector
  #The third element is the frequency of the third arrangement 
  #in the complex inversion
  aaf_pop_A = append(aaf_pop_A, 1-sum(aaf_pop_A))
  aaf_pop_B = append(aaf_pop_B, 1-sum(aaf_pop_B))
  
  #Estimate the expected heterozygosity
  expec_heterozygosity_A = 1 - sum(aaf_pop_A^2)
  expec_heterozygosity_B = 1 - sum(aaf_pop_B^2)
  
  #Estimate the heterozygosity in the sample (hs)
  h_s = (expec_heterozygosity_A + expec_heterozygosity_B) / 2
  
  #Estimate total heterozygosity (ht)
  ht_alleles = (aaf_pop_A + aaf_pop_B)/2
  #H = 1 - Σ(p^2)
  h_t = 1 - sum(ht_alleles^2)
  
  #Find the indexes of all ht equal 0
  h_t_0 = which(h_t == 0)
  
  #Fst = (Ht-Hs)/Ht
  fst = (h_t - h_s) / h_t
  print(h_t_0)
  #Replace those Nan values of fst with 0
  fst[h_t_0] = 0
  
  return(fst)
}


######################################################
# Compute Alternative Allele Frequencies
######################################################
compute_alt_af <- function (dataset){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  nNA = 0
  cN = c()
  alt_af = 0
  n_homozygous = 0
  n_heterozygous = 0
  n_NA = 0
  alt_af_array = c()
  for (i in c(1:nrows)){
    row = dataset[i,]
    n_homozygous = 0
    n_heterozygous = 0
    n_NA = 0
    #The columns are the individuals, so we count the number of copies of 
    #genotypes that have the allele '1'
    for (j in c(2:ncols)){
      genotype = toString(row[j])
      if(genotype == '0/1'){
        n_heterozygous = n_heterozygous + 1
      }else{
        if(genotype == '1/1'){
          n_homozygous = n_homozygous + 1
        }else{
          if(genotype == 'NA'){
            n_NA = n_NA + 1
          }
        }
      }
    }
    #Compute the frequency of the alternative allele ('1' allele)
    alt_af = ((n_homozygous * 2) + n_heterozygous) / ((ncols - 1 - n_NA) * 2)
    alt_af_array = c(alt_af_array, alt_af)
  }
  alt_af_array
}


######################################################
# Returns an array with the neutral and outlier
# data for a given population name
######################################################
get_data_by_population <- function(pop_name, p_neutral, p_outliers, p_lcols){
  l_cols_pop = grep(pop_name, p_lcols)
  l_cols_pop = c(1, l_cols_pop)
  pop_neutral = p_neutral[,l_cols_pop]
  pop_outliers = p_outliers[,l_cols_pop]
  return(list(pop_neutral,pop_outliers))
}

######################################################
#Remove rows (SNPs) with a population that does not 
#contain data for a minimum of individuals n_min_snails
######################################################
remove_incomplete_snps <- function(dataset, n_min_snails){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  n_samples = ncols - 1
  rows_keep = c()
  n_dataset = c()
  
  populations = unlist(strsplit(sample_names_analysis,"\\|"))
  col_numbers_populations = list()
  tmp_colnames_analysis = colnames(dataset)
  
  for(p in c(1:length(populations))){
    tmp_l_cols_pop = grep(populations[p], tmp_colnames_analysis)
    col_numbers_populations[[p]] = tmp_l_cols_pop
  }
  
  for (i in c(1:nrows)){
    row = dataset[i,]
    keep_SNP = TRUE
    #Examine if there are more than n_min_snails of NA snails 
    # within at least one population of the current time point
	#If that is the case, exclude the SNP
    for(p in c(1:length(populations))){
      tmp_l_cols_pop = unlist(col_numbers_populations[[p]])
      tmp_pop_row = row[tmp_l_cols_pop]
      #tmp_numdata_pop = tmp_pop_row[tmp_pop_row != 'NA']
      tmp_numdata_pop = tmp_pop_row[!is.na(tmp_pop_row)]
      if(length(tmp_numdata_pop) < n_min_snails){
        keep_SNP = FALSE
        print(paste("----> SNP ", i, "eliminated, # of Cols with data: ", length(tmp_numdata_pop), "/", length(tmp_l_cols_pop), "in POP ", populations[p]))
      }
    }
    
    if(keep_SNP){
      rows_keep = c(rows_keep, i)
    }
    
  }
  n_dataset = dataset[rows_keep, ]
  return(n_dataset)
}


##########################################################
#Remove columns (Snails) with percentage of of NAs
# greater than the indicated in the p_max_na_snps param
##########################################################
remove_incomplete_samples_by_percentage <- function(dataset, p_max_na_snps){
  nrows = dim(dataset)[1]
  ncols = dim(dataset)[2]
  cols_keep = c(1)
  n_dataset = c()
  #Samples must have info for > p_max_na_snps% of the SNPs
  threshold_na = trunc(nrows * p_max_na_snps) #Allow as much as p_max_na_snps% of NA SNPs per sample
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
#Remove rows (SNPs) with  >= p_max_na_snails % of snails with NA genotypes
##################################################################
remove_incomplete_snps_by_percentage <- function(dataset, p_max_na_snails){
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
    #Examine if there are more than p_max_na_snails% of NA snails for the SNP in total
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
    
  }
  n_dataset = dataset[rows_keep, ]
  return(n_dataset)
}


######################################################
# Translate the names of the snails in the shape data
# to the standard pop names in most analysis
# group_ref is a flag that indicates weather crab and wave should
# be grouped as single reference population
######################################################
translate_shape_names <- function(l_individuals, group_ref){
  tmp_translated_pop_name_list = c()
  for(ind in l_individuals){
    tmp_last_index = tail(unlist(gregexpr('_', ind)), n=1)
    tmp_pop_name = substr(ind,1,tmp_last_index-1)
    
    translated_pop_name = ""
    crab_18 = "Crab 2018"
    crab_21 = "Crab 2021"
    wave_18 = "Wave 2018"
    wave_21 = "Wave 2021"
    
    if(group_ref){
      crab_18 = "Crab"
      crab_21 = "Crab"
      wave_18 = "Wave"
      wave_21 = "Wave"
    }
    
    translated_pop_name = switch(  
      tmp_pop_name,  
      "Crab92"  = "Crab 1992",  
      "DON_92"  = "Crab 1992", 
      "Skerry96"= "Skerry 1996",
      "Skerry02"= "Skerry 2002",
      "Exp05"   = "Skerry 2005",
      "RamshC"  = crab_18,
      "Rskerry" = "Skerry 2018",
      "SkareW"  = wave_18,
      "DON_21"  = crab_21,
      "EXP_21"  = "Skerry 2021",
      "W_REF_21"= wave_21
    )  
    
    tmp_translated_pop_name_list = append(tmp_translated_pop_name_list, translated_pop_name)
  }
  
  return(tmp_translated_pop_name_list)
}


##################################################################
# Normalize the columns of a dataframe by minimum and maximum within
# the column
##################################################################
normalize_df_cols <- function(df_x) {
  norm_df_x = data.frame(matrix(nrow = nrow(df_x), ncol = 0))
  for(c in 1:ncol(df_x)){
    tmp_cx = df_x[,c]
    norm_cx = ((tmp_cx - min(tmp_cx, na.rm = TRUE)) / (max(tmp_cx, na.rm = TRUE) - min(tmp_cx, na.rm = TRUE)))
    norm_df_x = cbind(norm_df_x, norm_cx)
  }
  
  colnames(norm_df_x) <- colnames(df_x)
  return (norm_df_x)
}


##################################################################
# Translate a list of names of individuals (l_individuals) from
# a previous naming format into the standardized format for
# population names used in the analysis of the skerry
##################################################################
translate_shell_length_pop_names <- function(l_individuals, p_group_populations){
  tmp_translated_list = c()
  
  for(ind_name in l_individuals){
    
    tmp_last_index = tail(unlist(gregexpr('_', ind_name)), n=1)
    tmp_pop_name = substr(ind_name,1,tmp_last_index-1)
    
    tmp_translated = ""
    
    #If p_group_populations was set to true, the crab ecotype populations
    #will be grouped in Donor Crab, otherwise the year will be conserved
    if(p_group_populations){
      tmp_translated = switch(  
        tmp_pop_name,  
        "Crab92"  = "Donor Crab",  
        "DON_92"  = "Donor Crab", 
        "Skerry96"= "Skerry 1996",
        "Skerry02"= "Skerry 2002",
        "Exp05"   = "Skerry 2005",
        "RamshC"  = "Donor Crab",
        "Rskerry" = "Skerry 2018",
        "SkareW"  = "Reff. Wave",
        "DON_21"  = "Donor Crab",
        "EXP_21"  = "Skerry 2021",
        "W_REF_21"= "Reff. Wave"
      )
    }else{
      tmp_translated = switch(  
        tmp_pop_name,  
        "Crab92"  = "Skerry 1992",  
        "DON_92"  = "Skerry 1992", 
        "Skerry96"= "Skerry 1996",
        "Skerry02"= "Skerry 2002",
        "Exp05"   = "Skerry 2005",
        "RamshC"  = "Crab 2018",
        "Rskerry" = "Skerry 2018",
        "SkareW"  = "Wave Wave",
        "DON_21"  = "Crab 2021",
        "EXP_21"  = "Skerry 2021",
        "W_REF_21"= "Wave Wave"
      )  
    }
    
    tmp_translated_list = append(tmp_translated_list, tmp_translated)
  }
  return(tmp_translated_list)
}