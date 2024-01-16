################# Inversions Drift Simulations ###############################
# Based on the inversion frequencies in 1992 simulate the trajectories 
# over 30 years under a regime of drift and migration from the Wave ecotype.
# Plot the simmulated trajectories along with the actual data trajectory.
#
# Figure S12
#
#author: "Diego Garcia"
#date: "2023-11-03"
###############################################################################

#Set the working directory
setwd('.')

#Load libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)

#Fixed values
f = 2 #Number of generations per year based on the demographic inference
Y = 29 #Number of years of the experiment
n_repetitions = 100 #How many loci will be evolved by drift

#Column number in the 500K random draws file from the likelihood surface
c_N_0 = 1 #Starting population size (skerry) in ln scale
c_r = 2 #Growth rate in ln scale
c_K = 3 #Carrying capacity in ln scale
c_M = 4 #Migration rate (Haploid genomes, 1.6 snails) ln[M + 0.5]

#AF Trajectories of the inversions
#Frequencies of the inversions in Skerry in 1992, 2005, 2018 & 2021
iskerry_frequency_file <- "../Data/Inversions_Trajectories_Frequency_Skerry.txt"

#Frequencies of the inversions in combined samples (2018+2021) of the Crab and Wave ecotype populations
iwave_frequency_file <- "../Data/Inversions_Trajectories_Frequency_WaveCrab.txt"

#500K random draws of parameter combinations from the likelikhood surface
random_draws_LS_file  <- "../Data/Interpolation demographic parameters_random values 23 Sept 2023.csv"

#Load the datasets
iskerry_frequency_data <- read.table(iskerry_frequency_file, header = T, check.names = FALSE, sep = "\t")
iwave_frequency_data <- read.table(iwave_frequency_file, header = T, check.names = FALSE, sep = "\t")
random_draws_LS_data <- read.csv(random_draws_LS_file, header = F, check.names = FALSE)

#Reverse the log transformation of the parameters from the likelihood surface analysis
df_ls_parameters = exp(random_draws_LS_data[,c(1:3)])
df_ls_parameters = cbind(df_ls_parameters, exp(random_draws_LS_data[,c_M])-0.5)
colnames(df_ls_parameters) = c('N0','r','K','M')

list_inversions_plts = list() #Temporal list to add the plots for each inversion
plot_avg_drift = TRUE #When True adds a dashed line with the average drift

#Iterate over all inversion names
all_inversions = unique(iskerry_frequency_data$Inversion)
for(s_inv in all_inversions){
  
  #Subset the frequency trajectory of an inversion in the Skerry population
  s_inv_data <- iskerry_frequency_data[iskerry_frequency_data$Inversion == s_inv,]
  
  #Subset the frequency of an inversion in the Wave population
  w_inv_data <- iwave_frequency_data[iwave_frequency_data$Inversion == s_inv &
                                       iwave_frequency_data$Year == 'Wave',]
  
  #Identify the initial frequency of the Skerry (1992)
  tmp_inv_iaf <- s_inv_data[s_inv_data$Year == 'C 1992', 'Frequency']
  tmp_inv_colors <- s_inv_data[s_inv_data$Year == 'C 1992', 'Color']
  
  #The 'Color' column represent the arrangement number.
  #So, Complex inversions have more than one arrangement: Wave and Crab
  c_colors = c('Wave')
  if(length(tmp_inv_colors) > 1){
    c_colors = c('Wave', 'Crab')
  }
  
  #Create the drift plot of the inversion and add it to the list
  tmp_inv_plt <- prepare_drift_plot(tmp_inv_iaf, c_colors, s_inv, s_inv_data, w_inv_data)
  list_inversions_plts[[s_inv]] <- tmp_inv_plt
}

#Create a PDF File with the output of the plot
pdf('Inversions_Drift_Simulation_SupplementaryFigure_2024.pdf', width=7, height=8)
plot_grid(list_inversions_plts[['LGC1.1'   ]],
          list_inversions_plts[['LGC1.2'   ]],
          list_inversions_plts[['LGC2.1'   ]],
          list_inversions_plts[['LGC4.1'   ]],
          list_inversions_plts[['LGC6.1/2' ]],
          list_inversions_plts[['LGC7.1'   ]],
          list_inversions_plts[['LGC7.2'   ]],
          list_inversions_plts[['LGC9.1'   ]],
          list_inversions_plts[['LGC10.1'  ]],
          list_inversions_plts[['LGC10.2'  ]],
          list_inversions_plts[['LGC11.1'  ]],
          list_inversions_plts[['LGC14.1/2']],
          list_inversions_plts[['LGC17.1'  ]],
          nrow = 5,ncol = 3)

dev.off()



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
# Function that prepares the AF simulation under drift 
# of multiple replicates and generates the plot
####################################################################
prepare_drift_plot <- function(p_list_af_initial, p_list_colors, p_title, p_inv_data, p_w_inv_data){
  
  #Add the number of generations corresponding to each year of the trajectory
  p_inv_data$NewYear <- as.numeric(substr(p_inv_data$Year, 3, 6))
  tmp_inv_generations <- transform_year_to_generation(p_inv_data)
  p_inv_data$Generation <- tmp_inv_generations
  df_arrangement_trajectory <- p_inv_data[,c('NewYear', 'Generation','Frequency', 'Color')]
  
  #Modify the Color column to add the "Wave" or "Crab" color
  df_arrangement_trajectory$Color = rep(p_list_colors, each=length(unique(p_inv_data$Year)))
  
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
  df_wave_frequency$Allele <- 1
  df_wave_frequency$Rep <- 1
  df_wave_frequency$Generation <- max_generation + 10 #Just to add the wave after 2021
  
  # Randomly sample  WITHOUT replacement the sets of parameters to simulate
  # on each repetition (simulated trajectory):
  set.seed(5)
  random_rows_parameters = sample(1:nrow(df_ls_parameters), n_repetitions, replace=FALSE)
  
  df_drift_simulation = data.frame(matrix(nrow = 0, ncol = 3))
  #Complex inversions have two initial allele frequencies, so simulate two alleles
  tmp_allele = 100
  aux_i_allele = 0
  for(p_af_initial in p_list_af_initial){
    aux_i_allele = aux_i_allele + 1
    #Obtain the frequency of the arrangement in wave to simulate migration
    tmp_af_wave = iwave_frequency_data[iwave_frequency_data$Inversion == p_title, 'Frequency'][aux_i_allele]
    
    for(nr in seq(1:n_repetitions)){
  
      tmp_set_parameters = df_ls_parameters[random_rows_parameters[nr],]
      
      tmp_drift_af = drift_evolve_exclude_migrants(p_af_initial, tmp_af_wave, tmp_set_parameters)
      tmp_drift_af = cbind(tmp_drift_af, rep(nr,dim(tmp_drift_af)[1]), rep(tmp_allele,dim(tmp_drift_af)[1]), rep(p_list_colors[aux_i_allele],dim(tmp_drift_af)[1]))
      df_drift_simulation = rbind(df_drift_simulation, tmp_drift_af)
    }
    tmp_allele = tmp_allele + 1
  }
  
  colnames(df_drift_simulation) <- c("AF","Generation","Rep", "Allele", "Color")
  
  plt <- plot_lines(df_drift_simulation, df_arrangement_trajectory, df_wave_frequency, p_title)
  return(plt)
}




####################################################################
# Generate a lines plot with the AF of multiple repetitions
####################################################################
plot_lines <- function(p_drift_df, p_arrangement_trajectory_df, p_wave_frequency_df, p_title){
  palette_all <- c('Wave' = '#85C9E8', 'Crab' = '#E29949', '100' = '#d9ddde', '101' = '#d9ddde')
  
  
  data_drift <- p_drift_df
  drift_size_line = 0.08
  drift_alpha_line = 0.3
  drift_type_line = "solid"
  
  drift_size_line_avg = 0.5
  drift_alpha_line_avg = 1
  drift_type_line_avg = "dashed"
  
  x_labels <- p_arrangement_trajectory_df$NewYear
  x_labels <- append(x_labels, 'Wave')
  x_breaks <- p_arrangement_trajectory_df$Generation
  x_breaks <- append(x_breaks, max(p_wave_frequency_df$Generation))
  
  data_drift_avg <- p_drift_df %>%
    group_by(Allele, Color, Generation) %>%
    summarize(AF = mean(as.numeric(AF), na.rm = TRUE), Rep = min(Rep))
  
  plt_selected <- ggplot(data=data_drift, 
                         aes(x=as.numeric(Generation), 
                             y=as.numeric(AF), 
                             colour = factor(Color),
                             group = interaction(Rep,Allele)
                         )) + 
    geom_line(size = drift_size_line, alpha = drift_alpha_line, linetype = drift_type_line) + 
    #geom_line(aes(y=means), size=0.3, color="black")+
    #geom_point(size = 2)+
    theme_minimal()+
    #ylim(0,0.2)+
    theme(legend.position="none", 
          legend.title=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.3),
          plot.title = element_text(hjust = 0.5, size=8, face='bold'),
          axis.text=element_text(size=5),
          axis.title=element_text(size=5))+
    ggtitle(label = p_title)+
    ylim(0,1)+
    labs(x="", y = "Frequency")
  
  #Add the line of the average drift
  plt_selected <- plt_selected + 
    geom_line(data = data_drift_avg, 
              aes(x=as.numeric(Generation), 
                  y=AF, 
                  colour = factor(Color), 
                  group = factor(Color)),
              size = drift_size_line_avg, alpha = drift_alpha_line_avg, linetype = drift_type_line_avg) 
  
  #Add lines of the inversion trajectories
  plt_selected <- plt_selected + 
    geom_line(data = p_arrangement_trajectory_df, 
              aes(x=as.numeric(Generation), 
                  y=as.numeric(Frequency), 
                  colour = factor(Color), 
                  group = factor(Color)),
              size = 0.5) +
    geom_point(data = p_arrangement_trajectory_df, 
               aes(x=Generation, 
                   y=Frequency, 
                   colour = factor(Color), 
                   group = factor(Color)),
               size = 0.7)
  
  plt_selected <- plt_selected + 
    geom_point(data = p_wave_frequency_df, 
               aes(x=as.numeric(Generation), 
                   y=as.numeric(Frequency), 
                   colour = factor(Color)),
               size = 0.8) +
    scale_color_manual(values = palette_all)
  
  
  # custom X axis:
  plt_selected <- plt_selected + scale_x_continuous( label = x_labels, breaks= x_breaks ) 
  
  plt_selected
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