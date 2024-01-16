################################ Shell Length Trajectory #######################
#Generate the plots that show the trajectory of four phenotypic traits 
#in skerry and wave over time
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################

#Define the path and name of the files that contain phenotypic data
beige_file <- "../data/Beige.txt"
ridges_file <- "../data/Ridges.txt"
tesselation_file <- "../data/Tesselation.txt"
thickness_file <- "../data/Thickness.txt"

#Read and load the data files that contains phenotypic data
beige <- read.csv(beige_file, header = T, check.names = FALSE, na.strings="NA", row.names=NULL)
ridges <- read.csv(ridges_file, header = T, check.names = FALSE, na.strings="NA", row.names=NULL)
tesselation <- read.csv(tesselation_file, header = T, check.names = FALSE, na.strings="NA", row.names=NULL)
thickness <- read.csv(thickness_file, header = T, check.names = FALSE, na.strings="NA", row.names=NULL)

#Set the reference populations as only 'Wave' or 'Crab'
beige[beige$Sample == 'C 2021' | beige$Sample == 'C 2018',]$Sample <- 'Crab'
beige[beige$Sample == 'W 2021' | beige$Sample == 'W 2018',]$Sample <- 'Wave'
ridges[ridges$Sample == 'C 2021' | ridges$Sample == 'C 2018',]$Sample <- 'Crab'
ridges[ridges$Sample == 'W 2021' | ridges$Sample == 'W 2018',]$Sample <- 'Wave'
tesselation[tesselation$Sample == 'C 2021' | tesselation$Sample == 'C 2018',]$Sample <- 'Crab'
tesselation[tesselation$Sample == 'W 2021' | tesselation$Sample == 'W 2018',]$Sample <- 'Wave'
thickness[thickness$Sample == 'C 2021' | thickness$Sample == 'C 2018',]$Sample <- 'Crab'
thickness[thickness$Sample == 'W 2021' | thickness$Sample == 'W 2018',]$Sample <- 'Wave'

#Find the proportion of beige snails in the population
data_beige <- beige %>%
  group_by(Sample) %>%
  summarise(num_beige = sum(Beige, na.rm = TRUE), count = sum(Total, na.rm = TRUE)) %>%
  mutate(frequency = num_beige/count, type_pop = Sample)
  
#Find the proportion of ridged shells with respect to the donor population
data_ridges <- ridges %>%
  group_by(Sample) %>%
  summarise(num_ridged = sum(Ridges, na.rm = TRUE), count = n()) %>%
  mutate(frequency = num_ridged/count, type_pop = Sample)

#Find the proportion of UNtesselated shells in the population
data_tesselation <- tesselation %>%
  group_by(Sample) %>%
  summarise(num_tesselated = sum(Tesselation, na.rm = TRUE), count = n()) %>%
  mutate(frequency = 1-(num_tesselated/count), type_pop = Sample)

#Find the proportion of the change in thickness with respect to the donor population
expectation_thickness = mean(thickness[thickness$Sample == 'C 1992',]$Thickness, na.rm = TRUE)
data_thickness <- thickness %>%
  group_by(Sample) %>%
  summarise(mean_thickness = mean(Thickness, na.rm = TRUE),
            sd = sd(Thickness/expectation_thickness, na.rm = TRUE)) %>%
  mutate(frequency = mean_thickness/expectation_thickness, type_pop = Sample)

#Set the population name as 'Skerry' to all samples whose population is not 'Wave' nor 'Crab'
data_beige[data_beige$Sample != 'Crab' & data_beige$Sample != 'Wave', ]$type_pop <- 'Skerry'
data_ridges[data_ridges$Sample != 'Crab' & data_ridges$Sample != 'Wave', ]$type_pop <- 'Skerry'
data_tesselation[data_tesselation$Sample != 'Crab' & data_tesselation$Sample != 'Wave', ]$type_pop <- 'Skerry'
data_thickness[data_thickness$Sample != 'Crab' & data_thickness$Sample != 'Wave', ]$type_pop <- 'Skerry'

#Sort the phenotypes according to the populatio name
data_beige <- data_beige[order(data_beige$type_pop),]
data_ridges <- data_ridges[order(data_ridges$type_pop),]
data_tesselation <- data_tesselation[order(data_tesselation$type_pop),]
data_thickness <- data_thickness[order(data_thickness$type_pop),]

##########################################################
#Plot merged line plot of the phenotypes
#The output is a list of lists. 
# First there is a list with two elements: 
#   1-(plots of colour, digedgeness, and patterning)
#   2-(plots of thickness)
#Second, each sub element of the previous lists is a list with other two elements:
#   1-(Plots of Skerry) & 2-(Plots of Wave)
##########################################################
generate_phenotypes_plots <- function(){
  data_beige$trait <-  'Beige'
  data_ridges$trait <- 'Ridged shells'
  data_tesselation$trait <- 'Unpatterned shells'
  data_thickness$trait <- 'Thickness'
  
  #Merge the data of three phenotypic traits that will be plotted together
  data_merged <- rbind(data_beige[,c('Sample','frequency', 'type_pop','trait')],
                       data_ridges[,c('Sample','frequency', 'type_pop','trait')],
                       data_tesselation[,c('Sample','frequency', 'type_pop','trait')])
  
  data_merged_skerry <- data_merged[data_merged$Sample != 'Crab' & data_merged$Sample != 'Wave',]
  data_merged_skerry$year = as.numeric(substr(data_merged_skerry$Sample, 3,6))
  
  data_merged_wave <- data_merged[data_merged$Sample == 'Wave',]
  data_merged_wave$year = 'Wave'
  
  #Thickness is plotted separately
  data_thickness_skerry <- data_thickness[data_thickness$Sample != 'Crab' & data_thickness$Sample != 'Wave',]
  data_thickness_skerry$year = as.numeric(substr(data_thickness_skerry$Sample, 3,6))
  data_thickness_wave <- data_thickness[data_thickness$Sample == 'Wave',]
  data_thickness_wave$year = 'Wave'
  
  #Generate the plots of the first three traits for Skerry
  plt_phen_freq_skerry <- plot_merged_phenotype_trajectory(data_merged_skerry, '', FALSE, FALSE, 'Frequency')
  #Generate the plots of the first three traits for Wave
  plt_phen_freq_wave <- plot_merged_phenotype_trajectory(data_merged_wave, '', TRUE, FALSE, 'Frequency')
  plt_phen_freq_wave + guides(color = FALSE, size = FALSE)
  
  #Generate the plot of thickness for Skerry
  plt_phen_freq_tickness_skerry <- plot_merged_phenotype_trajectory(data_thickness_skerry, '', FALSE, TRUE, 'Relative to Crab')
  #Generate the plot of thickness for Wave
  plt_phen_freq_tickness_wave <- plot_merged_phenotype_trajectory(data_thickness_wave, '', TRUE, TRUE, 'Relative to Crab')
  
  #Arrange the output in a list, first the three traits
  list_phenotypes_frequencies_plt = list()
  list_phenotypes_frequencies_plt[['Skerry']] <- plt_phen_freq_skerry
  list_phenotypes_frequencies_plt[['Wave']] <- plt_phen_freq_wave
  
  #Arrange the output in a list, second the thickness
  list_phenotypes_thickness_plt = list()
  list_phenotypes_thickness_plt[['Skerry']] <- plt_phen_freq_tickness_skerry
  list_phenotypes_thickness_plt[['Wave']] <- plt_phen_freq_tickness_wave
  
  list_phenotypes = list(list_phenotypes_frequencies_plt, list_phenotypes_thickness_plt)
  
  return(list_phenotypes)
  
  #pdf('Phenotypic_Trajectory_Merged_down_up.pdf', width = 12, height = 3)
  #plot_grid(plt_skerry, plt_wave, nrow = 1, ncol = 2, rel_widths = c(0.9,0.1))
  #plot_grid(plt_tickness_skerry, plt_tickness_wave, nrow = 1, ncol = 2, rel_widths = c(0.9,0.1))
  #dev.off() 
  
}


############################################################################
# Creates a csv with the frequency or proportion of each phenotype for each year
############################################################################
write_frequencies_file <- function(){
  data_merged <- data.frame(Sample=data_beige$Sample)
  data_merged$f_Beige <- data_beige$frequency
  data_merged$f_Ridged <- data_ridges$frequency
  data_merged$f_Tesselated <- data_tesselation$frequency
  freq_thickenss <- append(data_thickness$frequency, NA, after = 1)
  freq_thickenss <- append(freq_thickenss, NA, after = 6)
  freq_thickenss <- append(freq_thickenss, NA, after = 8)
  data_merged$f_Thickness <- freq_thickenss
  write.csv(data_merged, 'FrequenciesPhenotypes.csv', quote = FALSE, row.names = FALSE)
}

######################################################
# Generates a line plot with frequency trajectories, all in one plot
######################################################
plot_merged_phenotype_trajectory <- function(data, p_title, p_hide_labels, p_error_bars, p_y_title){
  plt <- ggplot(data, 
               aes(x=if(p_hide_labels){factor(year)}else{year},
                   y=as.numeric(frequency), 
                   colour = trait, 
                   group = interaction(trait, type_pop)
                   )) + 
    geom_line(linewidth = line_width) + 
    geom_point(size = point_size, aes(shape = factor(trait)))+
    theme_minimal()+
    #xlim(1992,2021)+
    scale_color_manual(name = 'Trait', values = c('Beige' = "#D696CF", 
                                  'Ridged shells' = "#95CA81", 
                                  'Unpatterned shells' = "#E09157",
                                  'Thickness' = "#7A659E")) +
    scale_shape_manual(name = 'Trait', values = c('Beige' = 19, 
                                  'Ridged shells' = 15, 
                                  'Unpatterned shells' = 17,
                                  'Thickness' = 18)) +
    theme(legend.position=if(show_legend_geneal){"bottom"}else{'none'}, 
          legend.title=element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.3),
          plot.title = element_text(hjust = 0.5, size=12, face='bold'),
          axis.text=element_text(size=10),
          axis.title=element_text(size=10))+
    labs(title = p_title, x="", y = p_y_title)
  
  if(p_hide_labels){
    plt <- plt+theme(axis.text.y = element_blank())+
           labs(x="", y = element_blank())
            
  }else{
    plt <- plt+scale_x_continuous(breaks=c(1992, 1996, 2002, 2005, 2018, 2021))
  }
  
  if(p_error_bars){
    plt <- plt+geom_errorbar(aes(ymin=frequency-sd, ymax=frequency+sd), width=.2,
                  position=position_dodge(0.05))+
    ylim(0,1.2)
  }else{
    plt <- plt+ylim(0,1.2)
  }
  
  return(plt)
}

