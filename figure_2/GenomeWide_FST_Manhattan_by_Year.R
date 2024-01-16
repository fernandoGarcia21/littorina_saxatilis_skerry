################# Genome-wide FST Manhattan by Year ############################
#Generates a Manhattan plot for a specific year with the FST values of the 
#Skerry vs merged Crab ecotype population including both collinear loci and inversions. 
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################


#Define the path and name of the files that contain control and spatial outliers SNP genotypes
swedish_neutral_snp_file <- "../Data/SkerryExperiment_Neutral_NOLG12_Swedenfilter.txt"
swedish_outlier_snp_file <- "../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt"
#The status of the SNPs with respect to the expected range of the demographic inference
outside_envelope_outliers_file <- "../Data/SEQSNPTM006_OUTLIER_STATUS_OUT.txt"
outside_envelope_neutral_file <- "../Data/SEQSNPTM006_OUTLIER_STATUS_NEU.txt"

#Inversion frequency information
iSkerry_frequency_file <- "../Data/Inversions_Trajectories_Frequency_Skerry.txt"
iWaveCrab_frequency_file <- "../Data/Inversions_Trajectories_Frequency_WaveCrab.txt"
contig_pos_gm_file <- "../generals/contigPositionsGM.txt" #Coordinates of the contigs in the genetic map
inv_coordinates <- "../generals/Inv_Coordinates.txt" #Coordinates of the inversions (e.g. start and end in the genetic map)

#Read and load the data files with the genotypes
neutral <- read.table(swedish_neutral_snp_file, header = T, check.names = FALSE)
outliers <- read.table(swedish_outlier_snp_file, header = T, check.names = FALSE)
outside_envelope_outliers <- read.table(outside_envelope_outliers_file, header = T, check.names = FALSE)
outside_envelope_neutral  <- read.table(outside_envelope_neutral_file, header = T, check.names = FALSE)
iskerry_frequency_data <- read.table(iSkerry_frequency_file, header = T, check.names = FALSE, sep = "\t")
iWaveCrab_frequency_data <- read.table(iWaveCrab_frequency_file, header = T, check.names = FALSE, sep = "\t")
pos_gm_data <- read.table(contig_pos_gm_file, header = T, check.names = FALSE)
inv_coordinates_data <- read.table(inv_coordinates, header = T, check.names = FALSE)

#Define the keywords of the sample names (columns) that will be used in the analysis
sample_names_1992 = "DO.92"
sample_names_2005 = "SKR.05"
sample_names_2018 = "DonCrab|ExpSkerry|RefW"
sample_names_2021 = "DO.21|SRK.21|W.REF"

#Merge the column names that will be used in the analysis, this pattern
#is used to subset the columns
sample_names_analysis = paste(sample_names_1992, sample_names_2005, sample_names_2018,sample_names_2021,sep = "|")

#Only keep the SNPid column and the columns of actual genotypes (exclude extra info columns)
l_cols = colnames(neutral)
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


##############################################
# Creates a Manhattan plot of FST with collinear loci and inversions
# p_title is the title of the plot
# pop_skerry_col_labels is the pattern of the sample names for an specific year e.g. SKR-05
# pop_skerry_inversions_year is the Year in the skerry inversion frequencies file, e.g. S 2005
##############################################
generate_fst_manhattan_year <- function(p_title, pop_skerry_col_labels, pop_skerry_inversions_year){
  data_skerry = get_data_by_population(pop_skerry_col_labels, neutral, outliers, l_cols)
  
  data_crab = get_data_by_population("DonCrab|DO.21", neutral, outliers, l_cols)
  data_wave = get_data_by_population("RefW|W.REF", neutral, outliers, l_cols)
  
  ####################
  # Collinear
  ####################
  #CRAB AF
  crab_neutral_af = compute_alt_af(data_crab[[1]])
  crab_outlier_af = compute_alt_af(data_crab[[2]])
  
  #WAVE AF
  wave_neutral_af = compute_alt_af(data_wave[[1]])
  wave_outlier_af = compute_alt_af(data_wave[[2]])
  
  
  #Skerry AF
  skerry_neutral_af = compute_alt_af(data_skerry[[1]])
  skerry_outlier_af = compute_alt_af(data_skerry[[2]])
  
  ####################
  # Inversions
  ####################
  
  #Join the frequencies of the arrangements in crab to the frequencies of the arrangements in Skerry 
  tmp_Crab_frequency_data = iWaveCrab_frequency_data[iWaveCrab_frequency_data$Year == 'Crab',c('Inversion', 'Color', 'Frequency')]
  tmp_skerrycrab_frequency_data = left_join(iskerry_frequency_data, tmp_Crab_frequency_data, 
                                            by=c('Inversion', 'Color'),
                                            suffix=c('.Skerry','.Crab'))
  
  skerry_inversion_af = tmp_skerrycrab_frequency_data[tmp_skerrycrab_frequency_data$Year == pop_skerry_inversions_year, ]
  
  
  
  ######################## FST ######################################
  #Estimate the frequency of the wave allele in the skerry samples
  #And thus estimate the FST with respect to the average crab ecotype
  ####################################################################
  
  #Neutral
  fst_wave_allele_neutral = prepare_wave_af_and_FST_year(crab_neutral_af, skerry_neutral_af, wave_neutral_af, neutral, 'Control')
  
  #Outliers
  fst_wave_allele_outlier = prepare_wave_af_and_FST_year(crab_outlier_af, skerry_outlier_af, wave_outlier_af, outliers, 'Outlier')
  
  #Merge fst collinear datasets
  fst_wave_allele_collinear = rbind(fst_wave_allele_neutral, fst_wave_allele_outlier)
  
  
  #Inversions
  fst_wave_allele_inversion = prepare_inversion_fst_year(skerry_inversion_af)
  
  #Add a column only with the contig name and pos in the contig
  contig_snp = str_split_fixed(fst_wave_allele_collinear$SNP, "_", 3)[,c(1,2)]
  contig_snp =  apply( contig_snp[ , c(1,2) ] , 1 , paste , collapse = "_" )
  fst_wave_allele_collinear$contig_snp = contig_snp
  
  aux <<- fst_wave_allele_collinear
  
  ######################### GENERATE MANHATTAN PLOTS OF COLLINEAR AND INVERSIONS ######################
  plt_year = prepare_fst_plot_selected_loci(fst_wave_allele_collinear, fst_wave_allele_inversion, p_title)
  return(plt_year)
}


#############################################
# Prepares the data frame to generate a Manhattan plot.
# Classify the SNPs as selected or no selected according to 
# the expected range of the demographic inference.
# There are a few control loci that were identified as selected.
#############################################
prepare_fst_plot_selected_loci <- function(fst_wave_allele_data, p_inversions_fst, p_title){
  
  merged_envelope_data <- rbind(outside_envelope_outliers[,c('SNP','Status')], outside_envelope_neutral[,c('SNP','Status')])
  merged_envelope_data$Selected <- ifelse(merged_envelope_data$Status == TRUE, "Yes", "No")
  
  fst_wave_allele_data$selected = rep(NA, dim(fst_wave_allele_data)[1])
  fst_wave_allele_data$selected <- merged_envelope_data[match(fst_wave_allele_data$contig_snp, merged_envelope_data$SNP),'Selected']
  fst_wave_allele_data[is.na(fst_wave_allele_data$selected),]$selected <- 'No'
  
  #Order by chromosome and position
  fst_sorted = fst_wave_allele_data[
    order( fst_wave_allele_data[,'chr'], fst_wave_allele_data[,'pos'],  fst_wave_allele_data[,'pos_in_contig']),
  ]
  
  #Add the coordinates of the inversion segments to calculate the x axis
  tmp_df_cords_start_end = fst_sorted[,c('chr','pos')]
  tmp_df_cords_start_end = rbind(tmp_df_cords_start_end, 
                                 setNames(p_inversions_fst[,c('chr','inner_start')], names(tmp_df_cords_start_end)),
                                 setNames(p_inversions_fst[,c('chr','inner_end')], names(tmp_df_cords_start_end)))
  
  #Re-estimate the position on the LG of each SNP
  data_cum <- tmp_df_cords_start_end %>% 
    group_by(chr) %>% 
    summarise(max_bp = max(pos)) %>% 
    mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
    select(chr, bp_add)
  
  #Add the new position to the collinear FST dataset
  fst_sorted_new <- fst_sorted %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(pos_cum = pos + bp_add)
  
  
  #Add the new position to the inversions FST dataset
  inversions_fst_new <- p_inversions_fst %>% 
    inner_join(data_cum, by = "chr") %>% 
    mutate(inner_start_cum = inner_start + bp_add,
           inner_end_cum = inner_end + bp_add)
  
  # Prepare X axis
  #Add the cumulated coordinates of the inversion segments to calculate the x axis
  tmp_df_cords_start_end = fst_sorted_new[,c('chr','pos_cum')]
  tmp_df_cords_start_end = rbind(tmp_df_cords_start_end, 
                                 setNames(inversions_fst_new[,c('chr','inner_start_cum')], names(tmp_df_cords_start_end)),
                                 setNames(inversions_fst_new[,c('chr','inner_end_cum')], names(tmp_df_cords_start_end)))
  
  # Determine the coordinates of the shaded rectangles and the x position of the axis labels
  axis_fst_df <- tmp_df_cords_start_end %>% 
    group_by(chr) %>% 
    summarize(center=( max(pos_cum) + min(pos_cum) ) / 2 ,
              start=min(pos_cum),
              end=max(pos_cum))
  #Add the gray white intercalated
  axis_fst_df$color <- rep(c(1,0),8)
  
  ref_end_previouslg = append(0, axis_fst_df$end+0.0001)
  ref_end_previouslg = ref_end_previouslg[1:length(ref_end_previouslg)-1]
  axis_fst_df$start = ref_end_previouslg
  

  plt = plot_manhattan_selected(fst_sorted_new, axis_fst_df, inversions_fst_new, p_title)
  return(plt)
}


#############################################
# Creates a manhattan plot with as many colors
# as types of loci: eg. Neutral, Outliers, Inversions, True outliers
#############################################
plot_manhattan_selected <- function(fst_df, axis_fst_df, p_inversions_fst_segments, p_title){
  
  # Make the plot
  plt <- ggplot() +
    
    # Show all points
    geom_point(data=fst_df, aes(x=pos_cum, y=fst, color=as.factor(type), size = as.factor(selected), shape = as.factor(selected)), alpha=0.7) +
    #green
    #scale_color_manual(values = c(Inversion="#e07a5f", Neutral="#81b29a", Outlier="#3d405b")) +
    #blue
    scale_color_manual(values = c(Inversion="#1d3557", Control="#a8dadc", Outlier="#F34573",
                                  'Simple'="#1d3557", 'Complex'="#1d3557"), name=NULL) +
    scale_size_manual(values = c('Yes'=2, 'No'=1), name=NULL) +
    scale_shape_manual(name = "With evidence for selection:", values = c('Yes'=17, 'No'=1)) +
    # custom X axis:
    scale_x_continuous( label = axis_fst_df$chr, breaks= axis_fst_df$center , expand = c(0.01,0)) +
    #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Custom the theme:
    theme_minimal() +
    theme( 
      legend.position="none",
      legend.title = element_text(size = 8, face = "bold"),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_line(size = 0.1),
      plot.title = element_text(hjust = 0.5, face='bold', size=12,),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y=element_text(size=8),
      plot.margin = margin( 1,0.0,0.0,0.0, "cm")
      
    )+
    #ylim(0.0, 1.0)+
    labs(y = "FST", x="") +
    ggtitle(p_title)
  
  #Add the segments of the Inversions
  plt <- plt + geom_segment(data = p_inversions_fst_segments,
                            aes(x=inner_start_cum, xend=inner_end_cum, y = FST, yend = FST, 
                                color = as.factor(Type)), linewidth=0.8)
  
  #Gray rectangles to separate the chromosomes
  plt <- plt + geom_rect(data=axis_fst_df, mapping=aes(xmin=start, xmax=end, ymin=0, ymax=1, fill=as.factor(color)), alpha=0.15, show.legend = FALSE)+
    scale_fill_manual(values = c('1'="gray", '0'="white"), name=NULL)+
    guides(fill = FALSE)
  
  plt <- plt + geom_rect(data=p_inversions_fst_segments, mapping=aes(xmin=inner_start_cum, xmax=inner_end_cum, ymin=0, ymax=FST), fill = '#1D5B79', alpha=0.1, show.legend = FALSE)
  
  return(plt)
}


#############################################
# Creates a dataframe with the FST estimates for each inversion
# of the Skerry vs the Crab ecotype. The dataframe includes
# the coordinates of the inversion and the chromosome.
# Complex inversions are treated as 3-allelic loci (FST estimates)
#############################################
prepare_inversion_fst_year <- function(p_skerry_crab_inv_af){
  inversions_fst_df = data.frame(matrix(nrow = 0, ncol = 6))
  tmp_inversion_ids = unique(p_skerry_crab_inv_af$Inversion)
  tmp_is_simple_inversion = TRUE
  tmp_type_inversion = 'NA'
  
  #Create a dataframe with the coordinates of the inversions and FST values
  for(inv_id in tmp_inversion_ids){
    
    tmp_inv_data = p_skerry_crab_inv_af[p_skerry_crab_inv_af$Inversion == inv_id,]
    
    #If the inversion is simple
    if(dim(tmp_inv_data)[1] < 2){
      tmp_type_inversion = 'Simple'
      tmp_fst_skerry_crab = compute_fst(tmp_inv_data$Frequency.Skerry, tmp_inv_data$Frequency.Crab)
    }else{
      #Then the inversion is complex and the FST must be computed as a multiallelic locus
      tmp_is_simple_inversion = FALSE
      tmp_type_inversion = 'Complex'
      tmp_fst_skerry_crab = compute_fst_multiallelic(tmp_inv_data[,'Frequency.Skerry'], tmp_inv_data[,'Frequency.Crab'])
    }
    
    #Extract the coordinates of each inversion
    tmp_num_snps = max(tmp_inv_data$NumSnps)
    tmp_inversion_split = unlist(strsplit(unlist(strsplit(inv_id, 'LGC'))[2],'\\.'))
    tmp_lg = tmp_inversion_split[1]
    tmp_idx_inv = 1
    #If there are more than one inversion in a LG
    if(length(tmp_inversion_split) > 1){
      tmp_idx_inv = tmp_inversion_split[2]# "LGC14.1/2" = 1/2; LGC1.2 = 2
      tmp_idx_inv = unlist(strsplit(tmp_idx_inv, '/'))[1]
    }
    
    #This is one row DF with the start and end coordinates of the inversion
    tmp_coordinates_inv = inv_coordinates_data[inv_coordinates_data$LG == as.numeric(tmp_lg), ][as.numeric(tmp_idx_inv), ]
    
    inversions_fst_df = rbind(inversions_fst_df, data.frame('Inversion' = inv_id,
                                                            'Type' = tmp_type_inversion,
                                                            'FST' = tmp_fst_skerry_crab,
                                                            'chr' = tmp_coordinates_inv$LG,
                                                            'inner_start' = tmp_coordinates_inv$inner_start,
                                                            'inner_end' = tmp_coordinates_inv$inner_end,
                                                            'NumSnps' = tmp_num_snps))
    
  }
  
  return(inversions_fst_df)
}


#############################################
# CREATES A DATAFRAME WITH THE AF OF THE WAVE ALLELE
# IN A POPULATION SKERRY OR CRAB
# dataset_ref is used to extract the names of the contigs 
# and the positions in the genetic map
#############################################
prepare_wave_af_and_FST_year <- function(p_crab_af, p_skerry_af, p_wave_af, dataset_ref, type){
  
  num_snps = length(p_skerry_af)
  crab_af = p_crab_af[1:length(p_crab_af)]
  skerry_af = p_skerry_af[1:length(p_skerry_af)]
  wave_af = p_wave_af[1:length(p_wave_af)]
  
  #Plot only the most common allele in wave
  for(n in seq(1,num_snps)){
    #If the alternative allele in crab has a higher frequency
    #Than the alternative allele in wave, we must find the
    #frequency of the reference allele both in crab and skerry
    #this is the frequency of the most common allele in wave.
    if(crab_af[n] > wave_af[n]){
      skerry_af[n] = 1-skerry_af[n];
      crab_af[n] = 1-crab_af[n];
    }
  }
  
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
  
  ref_names = ref_names[with(ref_names, order(chr, pos)),]
  sorted_skerry = skerry_af[ref_names$row_number]
  sorted_crab = crab_af[ref_names$row_number]
  
  tmp_fst = compute_fst(sorted_skerry, sorted_crab)
  
  #create a data frame with the coordinate and af to plot
  data_af_year <- cbind(ref_names, sorted_skerry, sorted_crab, tmp_fst, rep(type,length(sorted_skerry)))
  colnames(data_af_year)[6:9] <- c('skerry', 'crab', 'fst', 'type')
  
  return(data_af_year)
}