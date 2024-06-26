################################ All Inversions Trajectories and Drift##########
#Generate the plots that show the trajectories of inversion arrangement frequencies
#of the skerry population from 1992 until 2021, as well as the combined Wave population
#
#author: "Diego Garcia"
#date: "2022-12-21"
################################################################################


#Fixed values
f = 2 #Number of generations per year
Y = 29 #Number of years of the experiment
n_repetitions = 1000 #How many loci will be evolved by drift

#Column number in the 5OOK random draws file from the likelihood surface
c_N_0 = 1 #Starting population size (skerry) in ln scale
c_r = 2 #Growth rate in ln scale
c_K = 3 #Carrying capacity in ln scale
c_M = 4 #Migration rate (Haploid genomes, 1.6 snails) ln[M + 0.5]

#Path and name of the input files
#These two files are generated by 1_InversionsGenerateTrajectoriesData.R
iskerry_frequency_file <- "../Data/Inversions_Trajectories_Frequency_Skerry.txt" #frequencies of the arrangements over time in Skerry
iwave_frequency_file <- "../Data/Inversions_Trajectories_Frequency_WaveCrab.txt" #frequencies of the arrangements over time in Wave and Crab

#500K random draws from the likelihood surface of the demographic inference
random_draws_LS_file  <- "../Data/Interpolation demographic parameters_random values 23 Sept 2023.csv" 

#Load data sets
iskerry_frequency_data <- read.table(iskerry_frequency_file, header = T, check.names = FALSE, sep = "\t")
iwave_frequency_data <- read.table(iwave_frequency_file, header = T, check.names = FALSE, sep = "\t")
iwave_frequency_data <- iwave_frequency_data[iwave_frequency_data$Year == 'Wave', ]
random_draws_LS_data <- read.csv(random_draws_LS_file, header = F, check.names = FALSE)

#Reverse the log transformation of the parameters from the likelihood surface analysis
df_ls_parameters = exp(random_draws_LS_data[,c(1:3)])
df_ls_parameters = cbind(df_ls_parameters, exp(random_draws_LS_data[,c_M])-0.5)
colnames(df_ls_parameters) = c('N0','r','K','M')


###########################################################
# Simulate drift for each simple and complex inversion,
# generate a arrangement frequency trajectory plot for each,
# and returns a list of plot objects
###########################################################
generate_inversions_plots <- function(){
  all_inversions = unique(iskerry_frequency_data$Inversion)
  
  simple_inversions = all_inversions[!all_inversions %in% c('LGC6.1/2','LGC14.1/2')]
  
  plt_simple_inversions = generate_plot_trajectories_group_inversions(simple_inversions, 'Simple Inversions', FALSE)
  plt_complex_inversion_6 = generate_plot_trajectories_group_inversions(c('LGC6.1/2'), 'Complex Inversion LGC6.1/2', TRUE)
  plt_complex_inversion_14 = generate_plot_trajectories_group_inversions(c('LGC14.1/2'), 'Complex Inversion LGC14.1/2', TRUE)
  
  list_inversions_plt = list(plt_simple_inversions, plt_complex_inversion_6, plt_complex_inversion_14)
  
  #pdf('Inversions/LS_SimpleInversionsAllInOne.pdf', width=10, height=6)
  #plot_grid(plt_simple_inversions, 
  #          plt_complex_inversion_6, 
  #          plt_complex_inversion_14,
  #          nrow = 3, ncol = 1, rel_heights = c(0.4,0.3,0.3))
  #dev.off()
  
  return(list_inversions_plt)
}




####################################################################
# Simulates drift for a group of inversions (e.g. all simple)
# and returns a ggplot object with the line plot of the trajectories
####################################################################
generate_plot_trajectories_group_inversions <- function(p_list_inversions, p_title, p_is_complex_inversion){
  
  all_inversions_trajectory_skr_df = data.frame(matrix(nrow = 0, ncol = 8))
  all_inversions_trajectory_wave_df = data.frame(matrix(nrow = 0, ncol = 6))
  all_quantiles_df = data.frame(matrix(nrow = 0, ncol = 4))
  tmp_last_gen = as.integer(f*Y)
  
  for(s_inv in p_list_inversions){
    print(paste('Inv:',s_inv))
    s_inv_data <- iskerry_frequency_data[iskerry_frequency_data$Inversion == s_inv,]
    w_inv_data <- iwave_frequency_data[iwave_frequency_data$Inversion == s_inv,]
    tmp_inv_iaf <- s_inv_data[s_inv_data$Year == 'C 1992', 'Frequency']
    tmp_inv_colors <- s_inv_data[s_inv_data$Year == 'C 1992', 'Color']
    
    inversion_trajectory_skr_df <- prepare_inv_skerry_df(s_inv_data)
    inversion_trajectory_skr_df$Inversion <- s_inv
    
    inversion_trajectory_wave_df <- prepare_inv_wave_df(inversion_trajectory_skr_df, w_inv_data)
    inversion_trajectory_wave_df$Inversion <- s_inv
    
    tmp_inv_quntiles <- determine_inversion_expectation(tmp_inv_iaf, s_inv, s_inv_data, w_inv_data)
    inversion_trajectory_skr_df$is_beyond_expectation <- FALSE
    inversion_trajectory_skr_df$Label <- 'Within neutral expectation'
    inversion_trajectory_skr_df$Allele <- 0
    
    indexes_last_generation = which(inversion_trajectory_skr_df$Generation == tmp_last_gen)
    idx_quantile = 1
    idx_fg = 1 #index of the first generation
    for(idx_lg in indexes_last_generation){
      tmp_last_year_af = inversion_trajectory_skr_df[idx_lg, 'Frequency']
      tmp_quant_arr = tmp_inv_quntiles[idx_quantile,]
      if(tmp_last_year_af < tmp_quant_arr$lowQuantil | tmp_last_year_af > tmp_quant_arr$upQuantil){
        inversion_trajectory_skr_df[seq(idx_fg,idx_lg),]$is_beyond_expectation <- TRUE
        
        #Set the label to show in the legend
        inversion_trajectory_skr_df[seq(idx_fg,idx_lg),]$Label <- inversion_trajectory_skr_df[seq(idx_fg,idx_lg),]$Inversion
        
        
        #Add a new column with an allele id (1 for first arrangement and 2 for second arrangement)
        inversion_trajectory_skr_df[seq(idx_fg,idx_lg),]$Allele <- idx_quantile
      }
      idx_quantile = idx_quantile + 1
      idx_fg = idx_lg+1
    }
    
    all_inversions_trajectory_skr_df = rbind(all_inversions_trajectory_skr_df, inversion_trajectory_skr_df)
    all_inversions_trajectory_wave_df = rbind(all_inversions_trajectory_wave_df, inversion_trajectory_wave_df)
    all_quantiles_df = rbind(all_quantiles_df,tmp_inv_quntiles)
  }
  
  #Add labels to show in the legend for the wave dataframe
  all_inversions_trajectory_wave_df$Label = unique(all_inversions_trajectory_skr_df[all_inversions_trajectory_skr_df$Inversion  %in% all_inversions_trajectory_wave_df$Inversion, c('Inversion', 'Allele', 'Label')])$Label
  all_inversions_trajectory_wave_df$NewYear = 'Wave'
  
  #Add a is_beyond_expectation column
  all_inversions_trajectory_wave_df$is_beyond_expectation = FALSE
  all_inversions_trajectory_wave_df[all_inversions_trajectory_wave_df$Label != 'Within neutral expectation',]$is_beyond_expectation = TRUE
  
  plt_inv_skerry <- plot_all_lines_trajectories(all_inversions_trajectory_skr_df,
                                                all_quantiles_df, p_title, FALSE, p_is_complex_inversion)
  plt_inv_wave<- plot_all_lines_trajectories(all_inversions_trajectory_wave_df,
                                             all_quantiles_df, p_title, TRUE, p_is_complex_inversion)
  
  list_inversion_plt = list()
  list_inversion_plt[['Skerry']] <- plt_inv_skerry
  list_inversion_plt[['Wave']] <- plt_inv_wave
  
  return(list_inversion_plt)
}


####################################################################
# Transforms the year data of the inversion trajectory into generations. 
####################################################################
transform_year_to_generation <- function(trajectory_data){
  list_year_generations = c()
  for (inv_y in trajectory_data$NewYear){
    #Convert the number of years into generations
    tmp_g = (inv_y - 1992) * f
    tmp_g = as.integer(tmp_g)
    list_year_generations = append(list_year_generations,tmp_g)
  }
  return(list_year_generations)
}



####################################################################
# Creates a dataframe with the trajectory of the inversions in skerry
####################################################################
prepare_inv_skerry_df <- function(p_inv_data){
  
  #Add the number of generations corresponding to each year of the trajectory
  p_inv_data$NewYear <- as.numeric(substr(p_inv_data$Year, 3, 6))
  tmp_inv_generations <- transform_year_to_generation(p_inv_data)
  p_inv_data$Generation <- tmp_inv_generations
  df_arrangement_trajectory <- p_inv_data[,c('NewYear', 'Generation','Frequency', 'Color')]
  
  return(df_arrangement_trajectory)
}


####################################################################
# Creates a dataframe with the frequency values of the inversions in Wave
####################################################################
prepare_inv_wave_df <- function(df_arrangement_trajectory, p_w_inv_data){
  
  #Create a DF for the wave data
  max_generation <- max(df_arrangement_trajectory$Generation)
  #For complex inversions choose the frequency and colors of two arrangements
  if(length(grep('LGC6|LGC14', p_w_inv_data)) > 0){
    df_wave_frequency <- p_w_inv_data
    df_wave_frequency$Color <- unique(df_arrangement_trajectory$Color)
  }else{
    df_wave_frequency <- p_w_inv_data[1,]
    df_wave_frequency$Color <- df_arrangement_trajectory[1,]$Color
  }
  
  df_wave_frequency$Allele <- seq(1,length(unique(df_arrangement_trajectory$Color)))
  df_wave_frequency$Rep <- 1
  df_wave_frequency$Generation <- max_generation + 5 #Just to add the wave after 2021
  
  return(df_wave_frequency)
}


####################################################################
# Function that prepares the AF simulation under drift 
# of multiple replicates and determines if the inversion
# frequency exceeds the 97.5% neutral envelope (quantiles)
####################################################################
determine_inversion_expectation <- function(p_list_af_initial, p_inv_name, p_inv_data, p_w_inv_data){
  
  #Complex inversions have two initial allele frequencies, so simulate two alleles
  aux_i_allele = 0
  tmp_last_gen = as.integer(f*Y)
  df_cols_percentiles = data.frame(matrix(nrow = 0, ncol = 4))
  
  # Randomly sample  WITHOUT replacement the sets of parameters to simulate
  # on each repetition (simulated trajectory):
  set.seed(5)
  random_rows_parameters = sample(1:nrow(df_ls_parameters), n_repetitions, replace=FALSE)
  
  for(p_af_initial in p_list_af_initial){
    aux_i_allele = aux_i_allele + 1
    #Obtain the frequency of the arrangement in wave to simulate migration
    tmp_af_wave = p_w_inv_data[,'Frequency'][aux_i_allele]
    
    #temporary save the simulated drift for each arrangement
    df_drift_simulation_arrangement = data.frame(matrix(nrow = 0, ncol = 3))
    for(nr in seq(1:n_repetitions)){
      
      tmp_set_parameters = df_ls_parameters[random_rows_parameters[nr],]
      
      tmp_drift_af = drift_evolve_exclude_migrants(p_af_initial, tmp_af_wave, tmp_set_parameters)
      tmp_drift_af = cbind(tmp_drift_af, rep(nr,dim(tmp_drift_af)[1]))
      df_drift_simulation_arrangement = rbind(df_drift_simulation_arrangement, tmp_drift_af)
    }
    
    colnames(df_drift_simulation_arrangement) <- c("AF","Generation","Rep")
    
    #Extract only the AF in the last generation of all replicates
    tmp_last_generation_afs <<- as.vector(df_drift_simulation_arrangement[df_drift_simulation_arrangement$Generation == tmp_last_gen, 'AF'])
    
    #obtain 2.5 and 97.5 percentiles of the distribution of the last generation simulations
    percentiles_interval = as.vector(quantile(tmp_last_generation_afs,na.rm = T,probs = c(0.025,0.975)))
    percentiles_interval = append(percentiles_interval, c(p_inv_name, aux_i_allele))
    
    df_cols_percentiles = rbind(df_cols_percentiles, percentiles_interval)
  }
  
  colnames(df_cols_percentiles) <- c('lowQuantil','upQuantil', 'Inversion', 'Allele')
  
  return(df_cols_percentiles)
  
}


####################################################################
# Generate a lines plot with the AF of all trajectories in the same plot
####################################################################
plot_all_lines_trajectories <- function(p_arrangement_trajectory_df, p_quantiles_df, p_title, p_hide_labels, p_is_complex){
  #palette_all <- c('1' = '#3EC70B', '2' = '#5463FF', '3' = '#EB4747', '4' = '#035397', '5' = '#E8630A', '6' = '#ED5EDD', '100' = '#d9ddde', '101' = '#d9ddde')
  print(p_arrangement_trajectory_df)
  #p_wave_frequency_df = p_wave_frequency_df[p_wave_frequency_df$Label != 'Within neutral expectation', ]
  
  palette_all <- c('TRUE' = '#FDA93F', 'FALSE' = '#686969')
  palette_labels <- c('LGC1.1'    = '#DC5964',
                      'LGC1.2'    = '#1E88E5',
                      'LGC2.1'    = '#E4DF22',
                      'LGC4.1'    = '#62E4E1',
                      'LGC6.1/2'  = '#C9AC90',
                      'LGC7.1'    = '#6782C3',
                      'LGC7.2'    = '#960083',
                      'LGC9.1'    = '#92D489',
                      'LGC10.1'   = '#DC5964',
                      'LGC10.2'   = '#1E88E5',
                      'LGC11.1'   = '#92D489',
                      'LGC14.1/2' = '#6782C3',
                      'LGC17.1'   = '#92D489',
                      'Within neutral expectation' = '#686969')
  
  
  trajectory_type_line = "solid"
  trajectory_alpha_line = 0.6
  
  tmp_years = c(1992,1996,2002,2005,2018,2021)
  tmp_generations = as.integer((tmp_years-1992)*f)
  
  x_labels <- p_arrangement_trajectory_df$NewYear
  #x_labels <- append(x_labels, 'Wave')
  x_breaks <- p_arrangement_trajectory_df$Generation
  #x_breaks <- append(x_breaks, max(p_wave_frequency_df$Generation))
  
  plt_selected <- ggplot(data=p_arrangement_trajectory_df, 
                         aes(x=if(p_hide_labels){factor(NewYear)}else{Generation},
                             y=Frequency, 
                             colour = if(p_is_complex){factor(Allele)}else{factor(Label)},
                             group = interaction(Inversion,Allele)
                         )) + 
    geom_line(aes(linewidth = factor(is_beyond_expectation),
                  linetype = factor(Allele)), 
              alpha = trajectory_alpha_line) + 
    geom_point(aes(size = factor(is_beyond_expectation)),
               shape=19,
               alpha = trajectory_alpha_line)+
    theme_minimal()+
    #ylim(0,0.2)+
    theme(legend.position=if(show_legend_geneal){"bottom"}else{'none'}, 
          legend.title=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.3),
          plot.title = element_text(hjust = 0.5, size=12, face='bold'),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10))+
    ggtitle(label = p_title)+
    ylim(0,1)+
    labs(x="", y = "Frequency")
  
  plt_selected <- plt_selected+
    scale_linewidth_manual(values = c('TRUE'=line_width, 'FALSE'=0.2))+
    scale_linetype_manual(values = c('0'='solid', '1'='dashed', '2'='solid'))+
    scale_size_manual(values = c('TRUE'=point_size, 'FALSE'=0.5))
  
  if(p_is_complex){
    #When the inversion is complex, the first arrangement is blue
    #because that is the color of the Wave arrangement
    plt_selected <- plt_selected+scale_color_manual(values = c('1'='#73BADA', '2'='#E09157'))
  }else{
    plt_selected <- plt_selected+scale_color_manual(values = palette_labels)
  }
  
  # custom X axis:
  if(p_hide_labels){
    plt_selected <- plt_selected+theme(axis.text.y = element_blank())+labs(x="", y = element_blank())
  }else{
    plt_selected <- plt_selected + scale_x_continuous( label = tmp_years, breaks= tmp_generations ) 
  }
  
  return(plt_selected)
}


####################################################################
# Function that models drift on a loci according to the initial AF
# And includes a bias by gene flow (M and AF in wave ecotype)
# but excludes the migrants from the population size
####################################################################
drift_evolve_exclude_migrants <- function(p_initial_af, p_wave_af, p_set_parameters){
  
  drift_af_df <- c(p_initial_af)
  
  total_generations = f * Y
  N_g = p_set_parameters$N0 #The size of the population, it will increase over each generation
  af_g = p_initial_af #The allele frequency at each generation 
  
  af_wave_migrants = p_wave_af * p_set_parameters$M #The frequency of the wave allele  is constant
  
  for(g in seq(1:total_generations)){
    
    #The population size at generation g with rate r
    N_g = population_logistic_growth(g, p_set_parameters)
    
    #Estimate the bias of the binomial distribution
    #as a proportion of the wave AF in the Skerry and the migrants
    expected_af_after_m = (af_g * (N_g - p_set_parameters$M) + af_wave_migrants) / N_g
    
    #Draw from a binominal distribution Ng haploid genomes with a bias of Skerry allele frequency
    tmp_genotypes_skerry_g = rbinom(n = N_g, size = 1, prob = expected_af_after_m)
    
    af_g = sum(tmp_genotypes_skerry_g) / N_g
    drift_af_df <- append(drift_af_df, af_g)
    
  }
  drift_af_df = cbind(drift_af_df, seq(0,total_generations))
  
  return(drift_af_df)
}

######################################################################
# Function that estimates the size of the population at generation p_T
# according to the logistic growth in the neutral model of the  
# demographic inference
######################################################################
population_logistic_growth <- function(p_T, p_set_parameters){
  n_T = 0
  tmp_divisor = 1+((p_set_parameters$K - p_set_parameters$N0)/p_set_parameters$N0)*exp(-p_set_parameters$r*p_T)
  n_T = p_set_parameters$K / tmp_divisor
  return (n_T)
}
