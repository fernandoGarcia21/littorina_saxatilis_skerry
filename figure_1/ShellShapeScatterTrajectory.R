################################ Shell Shape Scatter Trajectory ################
#Generate the scatter plots that show the trajectory of the shell shape 
#in skerry (including 1992) and wave over time
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################

#Read and load the data file that contains shape data
shape_phenotypes_data_file = '../Data/AllShapePhenotypicData.csv'
shape_phenotypes_data_original <- read.csv(shape_phenotypes_data_file, header = T, check.names = FALSE)
shape_phenotypes_data <- na.omit(shape_phenotypes_data_original) #exclude NA data

consider_min_max = TRUE #Adjust all scatter plots according to the minimum and maximum for a variable in all samples
use_log_scale = TRUE #Transform quantitative traits like shell_length into log scale

#Define the variables on the X and Y axis and the colour variable
x_variable = 'ln(gh)'
y_variable = 'shell_length'
cv = 'colour' #Qualitative variable

#Define two appearance settings
point_size_scatter = 0.5 #size of the points in the scatter plot
#The colours of the points based on a qualitative variable
colors_fixed=c('0' = '#1e2a54', 
               '1'='#e85159', 
               '3' = '#4ca169', 
               '?' = '#94d4c4',
               "multiple bands" = 'purple', 
               "patchy" = '#57755f',
               "band" = "#fac800",
               "band?" = "#eb7d00",
               "Beige" = "#d696cf",
               "Black" = "#1e2a54",
               "Black/Beige" = "#696900",
               "Yellow" = "#fac800",
               "Olive" = "#027818",
               "Orange" = "#eb7d00",
               'white' = "gray",
               'White' = "gray",
               "Brown" = "#a38348",
               "Eroded" = "#c584e8",
               na.value = "#5bb6eb"
)

##########################################################################
#Plot a scatter with the x-y axes by correlated quantitative traits
# And color the points by a categorical trait
##########################################################################
generate_phenotypes_scatter_plots <- function(){
  tmp_plot_i = 0
  data_to_plot = shape_phenotypes_data_original
  data_to_plot$population = translate_shape_names(data_to_plot$snailID, TRUE)
  
  vec_pops_order = c('1992','1996','2002','2005','2018','2021','Wave')
  shape_phenotypes_skerry <- data_to_plot[grepl(paste(vec_pops_order, collapse = '|'), data_to_plot$population),]
  
  
  skerry_populations = unique(shape_phenotypes_skerry$population)
  skerry_populations = skerry_populations[order(skerry_populations, decreasing = FALSE)]
  
  #identify the limits xy of each cualitative variable
  x_limits = c(min(shape_phenotypes_skerry[,x_variable], na.rm = TRUE),
               max(shape_phenotypes_skerry[,x_variable], na.rm = TRUE))
  
  y_limits = c(min(shape_phenotypes_skerry[,y_variable], na.rm = TRUE),
               max(shape_phenotypes_skerry[,y_variable], na.rm = TRUE))
  
  if(use_log_scale && !x_variable %in% c('ln(gw)','ln(gh)')){
    x_limits = log(x_limits)
  }
  
  if(use_log_scale && !y_variable %in% c('ln(gw)','ln(gh)')){
    y_limits = log(y_limits)
  }
  
  tmp_list_years_plts <- list()
  is_first = TRUE
  #One column in a row is a population year (1996, 2002, 2005...)
  for(sp in skerry_populations){
    tmp_plot_i = tmp_plot_i + 1
    tmp_year_name = vec_pops_order[tmp_plot_i]
    #Subset the data for one population
    tmp_data_plot = shape_phenotypes_skerry[grepl(sp, shape_phenotypes_skerry$population),c(x_variable,y_variable,cv, 'snailID')]
    colnames(tmp_data_plot) <- c('x','y','col','ind')
    tmp_data_plot = na.omit(tmp_data_plot)
    
    
    if(use_log_scale && !x_variable %in% c('ln(gw)','ln(gh)')){
      tmp_data_plot$x = log(tmp_data_plot$x)
    }
    if(use_log_scale && !y_variable %in% c('ln(gw)','ln(gh)')){
      tmp_data_plot$y = log(tmp_data_plot$y)
    }
    
    #Create the plot
    plt <- create_plot_pair(tmp_data_plot,
                            sp,
                            x_variable, 
                            y_variable, 
                            x_limits,
                            y_limits,
                            cv,
                            sp,
                            is_first,
                            consider_min_max)
    
    tmp_list_years_plts[[tmp_year_name]] <- plt
    is_first = FALSE
  }
  
  return(tmp_list_years_plts)
}




###########################################################################
# Creates a ggplot object with the scatter plot of two variables (x_variable and y_variable)
###########################################################################
create_plot_pair <- function(p_data_plot, skerry_pop, 
                             x_variable, y_variable, 
                             x_limits, y_limits,
                             col_variable, p_title,
                             is_first, consider_min_max){

  plt <- ggplot(data = p_data_plot, aes(x=x, y=y, 
                                        color=factor(col),
                                        label=ind))+
    geom_point(size=point_size_scatter, alpha=0.6, shape=19)+
    theme_minimal()+
    theme(legend.title = element_blank(),
          legend.position=if(show_legend_geneal){"bottom"}else{'none'}, 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.3),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.1)
          )+ 
    scale_color_manual(values=colors_fixed)#+
    #ggtitle(p_title)
  
  #If the x and y axis should be adjusted for all plots to have the same scale
  if(consider_min_max){
    plt <- plt+xlim(x_limits)+ylim(y_limits)
  }
  
  #Only show x and y labels for the first scatter plot (1992)
  if(is_first){
    plt <- plt + xlab(x_variable) + ylab(y_variable)
  }else{
    plt <- plt + xlab(' ') + ylab(' ')#+
      #theme(axis.text=element_text(colour="white"))
      #scale_y_continuous(labels = function(breaks_y) {rep_along(breaks, "")})+
      #scale_x_continuous(labels = function(breaks_x) {rep_along(breaks, "")})
  }

  return(plt)
}
