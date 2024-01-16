################################ Shell Length Trajectory #######################
#Generate the plots that show the trajectory of the shell length 
#in skerry and wave over time
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################


#Read and load the data file that contains shape data
shape_data_file = '../Data/AllShapePhenotypicData.csv'
shape_data <- read.csv(shape_data_file, header = T, check.names = FALSE)

#Get the list of individuals in the shape data file
individuals_list = shape_data$snailID 

#Add a column with the name of the population after translating 
#the individuals names into a standardized format
shape_data$population = translate_shell_length_pop_names(individuals_list, FALSE)

# Split the name of the population of the individuals to obtain the year 
list_names_split = strsplit(shape_data$population, split = " ")
#The matrix is an auxiliary structure to easy manipulate the output of the strsplit
matrix_years = matrix(unlist(list_names_split),ncol=2,byrow=T)

#Create a dataframe with pop name, sell length, and year
df_shell_length <- data.frame(population = matrix_years[,1], 
                              length = shape_data$shell_length,
                              year = matrix_years[,2])

######################################################
# Generates a list with two plots: Skerry and Wave
######################################################
generate_shell_length_plots <- function(){
  #Estimate average shell length by year and its SD
  df_summary <- df_shell_length %>%
    group_by(population, year) %>%
    summarise(mean = mean(length, na.rm=TRUE),
              sd = sd(length, na.rm=TRUE))
  
  #Exclude Crab 2018 and Crab 2021
  df_shell_length_skerry <- df_summary %>%  filter(population =='Skerry')
  df_shell_length_wave <- df_summary %>%  filter(population =='Wave')
  
  plt_length_skerry <- plot_shell_length(df_shell_length_skerry, FALSE)
  plt_length_wave <- plot_shell_length(df_shell_length_wave, TRUE)
  
  #pdf('ShellLengthTrajectory.pdf', width=12, height=2.5)
  #plot_grid(plt_length_skerry, plt_length_wave, nrow = 1, ncol = 2, rel_widths = c(0.9,0.1))
  #dev.off()
  
  list_shell_length_plt = list()
  list_shell_length_plt[['Skerry']] <- plt_length_skerry
  list_shell_length_plt[['Wave']] <- plt_length_wave
  
  return(list_shell_length_plt)
}


######################################################
# Generate a scatter plot with the length of the shells
# and the limits of the SD
######################################################
plot_shell_length <- function(df_data, p_hide_labels){
  
  plt <- ggplot(df_data, aes(x=if(p_hide_labels){factor(year)}else{as.numeric(year)}, y=mean)) + 
    geom_point(size = point_size, colour='#FAC800', shape = 19)+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.5,
                  position=position_dodge(0.05), colour='#FAC800')+
    theme_minimal()+
    xlab("")+
    ylab("mm")+
    theme(
      legend.position=if(show_legend_geneal){"bottom"}else{'none'}, 
      legend.title=element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3),
      plot.title = element_text(hjust = 0.5, size=12, face='bold'),
      axis.text=element_text(size=10),
      axis.title=element_text(size=10))
  
  if(p_hide_labels){
    plt <- plt+theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
      labs(x="", y = element_blank())
    
  }else{
    plt <- plt+scale_x_continuous(breaks=c(1992, 1996, 2002, 2005, 2018, 2021))
  }
  
  plt <- plt+ylim(3.0,12)
  
  return(plt)
}
