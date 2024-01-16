################# Phenotypes Shape Correlations ###############################
# Analyse the quantitative and qualitative traits to find uncorrelated pairs 
# of traits and plot the scatter and barplots of the traits.
# This analysis permited to explore possible bimodalities along the time.
#
# Figure S3: Correlations and Box Plots
# Figure S4: Scatter Plots
# Figure S5: Barplots
#
#author: "Diego Garcia"
#date: "2023-06-10"
###############################################################################

#Set the working directory
setwd('.')

#Load libraries
library(matlib)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(dplyr)
library('NbClust')
source('../generals/GeneralsSkerry.R')

#Define the path and name of the input phenotypic file
shape_phenotypes_data_file = '../Data/AllShapePhenotypicData.csv'

#Load the phenotypes dataset from the file
shape_phenotypes_data_original <- read.csv(shape_phenotypes_data_file, header = T, check.names = FALSE)
shape_phenotypes_data <- na.omit(shape_phenotypes_data_original) #Exclude rows with NA values

#Keep only the necesary columns from the dataset (traits and population name)
quantitative_data <- shape_phenotypes_data[,c('ln(gw)', 'ln(gh)', 
                                              'r0', 'h0', 'a0', 
                                              'c', 'shell_length', 'avg_thickness',
                                              'population'
                                              )]
#extract the population names
tmp_pop_names <- quantitative_data$population

#Define the colour palet tho be used in the scatter and bar plots
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
#Generate a correlation plot among all phenotypic traits (Figure S3 A)
##########################################################################

#Subset the data by population
quantitative_data_crab <- quantitative_data[grepl('Crab', tmp_pop_names),]
quantitative_data_wave <- quantitative_data[grepl('Wave', tmp_pop_names),]
quantitative_data_skerry <- quantitative_data[grepl('Skerry', tmp_pop_names),]
col_corr = c(c("#73BADA", "white", "#E09157")) #Define the colour palette of the correlations

#Create a PDF output file
pdf('Shape_Correlations_SupplementaryFigure_2024.pdf', width = 10, height =4)
par(mfrow = c(1, 3)) #The output is a matrix of one row and three columns

#Estimate and plot the correlation between quantitative traits of all samples
corrplot(cor(quantitative_data[,1:8]), 
         method="circle", type="upper", 
         title = 'All Data (Crab, Wave, & Skerry)', 
         mar = c(2, 0, 1, 0), 
         colorRampPalette(col_corr)(100),
         tl.col = 'black'
         )
#Estimate and plot the correlation between quantitative traits of only Crab samples
corrplot(cor(quantitative_data_crab[,1:8]), 
         method="circle", type="upper", 
         title = 'Crab (1992, 2018, & 2021)', 
         mar = c(2, 0, 1, 0),
         colorRampPalette(col_corr)(100),
         tl.col = 'black')
#Estimate and plot the correlation between quantitative traits of only Wave samples
corrplot(cor(quantitative_data_wave[,1:8]), 
         method="circle", type="upper", 
         title = 'Wave (2018 & 2021)', 
         mar = c(2, 0, 1, 0),
         colorRampPalette(col_corr)(100),
         tl.col = 'black')
#corrplot(cor(quantitative_data_skerry[,1:8]), method="circle", type="upper", title = 'Skerry',mar = c(2, 0, 1, 0))
dev.off()
##########################################################################



##########################################################################
#Plot boxplots of all quantitative traits in Crab and Wave (Figure S3 B)
#########################################################################

data_to_plot = shape_phenotypes_data_original
#Add the name of the population according to sample name.
#Group the Crab and Wave samples (do not distinguish by year).
data_to_plot$population = translate_shape_names(data_to_plot$snailID, TRUE)

#Subset the data to keep only the columns (traits) that we want to plot
data_to_plot_quant <- data_to_plot[,c('ln(gw)', 'ln(gh)', 
                                              'r0', 'h0', 'a0', 
                                              'c', 'shell_length', 'avg_thickness',
                                              'population')]

#Subset the data to plot only Crab and Wave populations data, exclude Skerry
data_q_cw <- data_to_plot_quant[data_to_plot_quant$population == 'Crab' | data_to_plot_quant$population == 'Wave',]

#Create a dataframe with three columns: variable, value (y), population (group)
df_to_plot = create_xy_dataframe(data_q_cw)

#Create a PDF file with the output figure
pdf('ShapeBoxPlots_SupplementaryFigure_2024.pdf', width = 10, height = 6)

#Create a boxplot with Crab and Wave phenotypes
ggplot(data = df_to_plot, aes(x = variable, y = y, 
                              color=factor(group),
                              fill = factor(group)))+
  geom_boxplot(alpha = 0.1) +
  geom_point(position = position_dodge(width=0.75),
             alpha = 0.5)+
  facet_wrap(~variable, scale="free")+
  theme_minimal()+
  theme(legend.title = element_blank(),
        plot.title = element_blank(),
        panel.grid.major = element_blank()
        )+
  scale_color_manual(values=c('Wave' = '#73BADA', 'Crab' = '#E09157'))+
  scale_x_discrete()
dev.off()
##########################################################################



##########################################################################
#Plot a scatter with the x-y axes by uncorrelated quantitative traits
# And color the points by a categorical trait.
# Uncorrelated traits (e.g. ln(gh) vs thickness) were obtained by 
# visual validation of the corr plots and box plots.
# (Figure S4)
##########################################################################
data_to_plot = shape_phenotypes_data_original
data_to_plot$population = translate_shape_names(data_to_plot$snailID, TRUE)

#Subset the data to keep only the samples that we need to visualize. 
# E.g. exclude Skerry 2018 and 2021 because it was very similar to Wave already.
shape_phenotypes_skerry <- data_to_plot[grepl('Crab 1992|1996|2002|2005|Wave', data_to_plot$population),]

#Sort the population names ascending
skerry_populations = unique(shape_phenotypes_skerry$population)
skerry_populations = skerry_populations[order(skerry_populations, decreasing = FALSE)]

#Define the column names of the qualitative traits to colour the samples
cualitative_variables = c('colour','tesselated','ridged')

#Define the pairs of uncorrelated quantitative traits (x vs y)
cualitative_variables_correlated = data.frame(x=c('ln(gh)','ln(gh)',
                                                  'ln(gw)','ln(gw)'),
                                              y=c('avg_thickness','shell_length',
                                                  'avg_thickness','shell_length'))

num_pairs = dim(cualitative_variables_correlated)[1]
n_rows = length(cualitative_variables)
n_cols = length(skerry_populations)+1 #Add one col because of the labels (1st col)
#Estimate the width of the columns on a page
w_cols = 0.95/(n_cols-1)
w_list = c(0.05,rep(w_cols,(n_cols-1)))
consider_min_max = TRUE #If True, adjust the X and Y axis to be the same for all plots
use_log_scale = TRUE #If True, log transform the shell thickness and length.

#Create a PDF file with the output plot 
pdf('ScatterShapePhenotypes_reduced_skerry_Final_Supplement.pdf', width = 15, height =5)
#One page in the document is a pair of shape parameters
for(i_p in 1:num_pairs){
  #Arrange a grid per page
  x_variable = cualitative_variables_correlated[i_p,'x']
  y_variable = cualitative_variables_correlated[i_p,'y']
  tmp_list_cualitatives = list()
  tmp_list_years_plts <- list()
  tmp_plot_i = 0
  is_first = TRUE
  
  #identify the limits xy of each cualitative variable
  x_limits = c(min(shape_phenotypes_skerry[,x_variable], na.rm = TRUE),
               max(shape_phenotypes_skerry[,x_variable], na.rm = TRUE))
  
  y_limits = c(min(shape_phenotypes_skerry[,y_variable], na.rm = TRUE),
               max(shape_phenotypes_skerry[,y_variable], na.rm = TRUE))
  
  if(use_log_scale){
    y_limits = log(y_limits) #Log transform thickness and length
  }
  
  #One row in one page is a qualitative variable (color, tessellated, ridged)
  for(cv in cualitative_variables){
    tmp_plot_i = tmp_plot_i + 1
    tmp_list_years_plts[[paste(tmp_plot_i)]] <- generate_plot_text(cv)
    #One column in a row is a population year (1996, 2002, 2005...)
    for(sp in skerry_populations){
      
      #Subset the data for one population
      tmp_data_plot = shape_phenotypes_skerry[grepl(sp, shape_phenotypes_skerry$population),c(x_variable,y_variable,cv, 'snailID')]
      tmp_data_plot$shape = 'normal'
      colnames(tmp_data_plot) <- c('x','y','col','ind','shape')
      tmp_data_plot = na.omit(tmp_data_plot)
      
      if(use_log_scale){
        tmp_data_plot$y = log(tmp_data_plot$y)
      }

      #Create a ggplot object (scatter) with the x and y pair of variables
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
      
      tmp_plot_i = tmp_plot_i + 1
      tmp_list_years_plts[[tmp_plot_i]] <- plt
      
    }
    is_first = FALSE
  }
  
  #Create a text object of the variable pairs to plot as a page title
  title_page=text_grob(paste(x_variable, 'vs', y_variable), 
                       size = 10, face = "bold", lineheight=5)
  
  #Plot a page in the file
  grid.arrange(grobs = tmp_list_years_plts, 
               ncol = n_cols, nrow=n_rows, 
               widths=w_list,
               bottom = text_grob(paste('Page', i_p, '/', num_pairs)),
               top = title_page)
  
}
dev.off()




##########################################################################
#Plot a stacked barplot of binned quantitative traits 
# including colours according to % of qualitative traits
# Uncorrelated traits (e.g. ln(gh) vs thickness) were obtained by 
# visual validation of the corr plots and box plots.
# (Figure S5)
##########################################################################

data_to_plot = shape_phenotypes_data_original
data_to_plot$population = translate_shape_names(data_to_plot$snailID, TRUE)

#Subset the data to keep only the samples that we need to visualize. 
# E.g. exclude Skerry 2018 and 2021 because it was very similar to Wave already.
shape_phenotypes_skerry <- data_to_plot[grepl('1992|1996|2002|2005|Wave', data_to_plot$population),]

#Extract and sort the population names ascending
skerry_populations = unique(shape_phenotypes_skerry$population)
skerry_populations = skerry_populations[order(skerry_populations, decreasing = FALSE)]

#Define the column names of the qualitative traits to colour the samples
cualitative_variables = c('ridged','colour','tesselated')

#Define the pairs of uncorrelated quantitative traits (x vs y)
x_variables = c('shell_length','avg_thickness','ln(gw)','ln(gh)')

#Identify the number of rows each page
n_rows = length(cualitative_variables)
n_cols = length(skerry_populations)+1 #Add one col because of the labels (1st col)
#Estimate the width of the columns on a page
w_cols = 0.95/(n_cols-1)
w_list = c(0.05,rep(w_cols,(n_cols-1)))
use_log_scale = FALSE #If True, log transform shell thickness and length
show_frequency = FALSE #If True, show the frequencies on the bars
n_breaks_x_axis = 4 #Number of bins (bars) on the X axis
tmp_pdf_name = ''

if(use_log_scale){
  tmp_pdf_name = 'BarPlotShapePhenotypes_LOG_Final_Supplement'
}else{
  tmp_pdf_name = 'BarPlotShapePhenotypes_Final_Supplement'
}
pdf(paste0(tmp_pdf_name,'_',n_breaks_x_axis,'.pdf'), width = 15, height =8)

#One page in the document is a pair of shape parameters
for(i_p in 1:length(x_variables)){
  #Arrange a grid per page
  x_variable = x_variables[i_p]
  tmp_list_cualitatives = list()
  tmp_list_years_plts <- list()
  tmp_plot_i = 0
  is_first = TRUE
  
  #One row in one page is a cualitative variable (color, tessellated, ridged)
  for(cv in cualitative_variables){
    tmp_plot_i = tmp_plot_i + 1
    tmp_list_years_plts[[paste(tmp_plot_i)]] <- generate_plot_text(cv)
    
    #Summarize the counts of the variables to feed the barplot
    summary_ref_data = na.omit(shape_phenotypes_skerry[, c('population', x_variable, cv)])
    colnames(summary_ref_data) <- c('pop','x','col')
    if(use_log_scale && !x_variable %in% c('ln(gw)','ln(gh)')){
      summary_ref_data$x = log(summary_ref_data$x) #Log transform shell thickness and length
    }
    
    #Create a dataframe with the beans according the data summary
    grouped_summary_df = summary_ref_data %>% 
      mutate(x= cut(x, breaks=n_breaks_x_axis)) %>%
      group_by(pop, x, col) %>% 
      summarise(count = n())
    
    x_axis = levels(grouped_summary_df$x)
    x_axis_ref = data.frame(x=x_axis, id_x=seq(1:n_breaks_x_axis))
    
    
    #One column in a row is a population year (1996, 2002, 2005...)
    for(sp in skerry_populations){
      
      #subset the summary data for the actual population to plot
      grouped_data_plot = grouped_summary_df[grouped_summary_df$pop == sp, ]
      
      #Transform the bin names to continuos values (1,2,3...)
      tmp_id_x_list = c()
      for(tmp_x in grouped_data_plot$x){
        tmp_id_x_list = append(tmp_id_x_list, x_axis_ref[x_axis_ref$x == tmp_x,]$id_x)
      }
      grouped_data_plot$id_x = tmp_id_x_list
      
      #Create a barplot of a quantitative trait and coloured bars by a qualitatitive trait
      plt <- create_bar_plot_pair(grouped_data_plot,
                              sp,
                              x_variable, 
                              show_frequency,
                              x_axis,
                              cv,
                              sp,
                              is_first)
      
      tmp_plot_i = tmp_plot_i + 1
      tmp_list_years_plts[[tmp_plot_i]] <- plt
    }
    is_first = FALSE
  }
  
  #Create a title for the page
  title_page=text_grob(paste(x_variable), x = 0, hjust = 0,
                       size = 20, face = "bold", lineheight=5)
  
  #Plot the barplots on the page
  grid.arrange(grobs = tmp_list_years_plts, 
               ncol = n_cols, nrow=n_rows, 
               widths=w_list,
               bottom = text_grob(paste('Page', i_p, '/', length(x_variables))),
               top = title_page)
}
dev.off()


###############################################################
# Creates a scatter plot object from a dataframe
###############################################################
create_plot_pair <- function(p_data_plot, skerry_pop, 
                             x_variable, y_variable, 
                             x_limits, y_limits,
                             col_variable, p_title,
                             is_first, consider_min_max){
  text_title = ""
  if(is_first){
    text_title = p_title
  }
  
  plt <- ggplot(data = p_data_plot, aes(x=x, y=y, 
                                        color=factor(col),
                                        label=ind))+
    geom_point(size=2, alpha=0.6, shape=19)+
    #geom_text(hjust=0, vjust=0, size=2)+
    theme_minimal()+
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5))+ 
          #panel.grid.minor = element_blank())+
    scale_color_manual(values=colors_fixed)+
    scale_shape_manual(values=c(normal=19,highlight=17))+
    xlab(x_variable)+
    ylab(y_variable)+
    ggtitle(text_title)
  
  #If the x and y axis should be adjusted for all plots to have the same scale
  if(consider_min_max){
    plt <- plt+xlim(x_limits)+ylim(y_limits)
  }
  
  return(plt)
}

###################################################
# Creates a ggplot object with a text
###################################################
generate_plot_text <- function(message){
  plt <- ggplot() + annotate("text",
                    x = 0,
                    y = 0,
                    size = 5,
                    angle=90, 
                    fontface="bold",
                    label = message) + 
                    theme_void()
  return(plt)
}


#############################################################
#Create a dataframe with three columns:
#1: The variable name (x)
#2: The values of the variable (y)
#3: The group of the sample (e.g. the population)
#############################################################
create_xy_dataframe <- function(df_data){
  df_to_plot = data.frame(matrix(nrow = 0, ncol = 3,
                                 dimnames=list(NULL, c("variable", "y", "group"))))
  n_rows = nrow(df_data)
  for(col in 1:(ncol(df_data)-1)){
    tmp_df_column = data.frame(
      variable = rep(colnames(df_data)[col],n_rows),
      y = df_data[,col],
      group = df_data$population
    )
    df_to_plot = rbind(df_to_plot, tmp_df_column)
  }
  df_to_plot = na.omit(df_to_plot)
}


########################################################################
# Creates a barplot plot object for a dataframe
########################################################################
create_bar_plot_pair <- function(p_data_plot, skerry_pop, 
                             x_variable, 
                             show_freq,
                             x_axis,
                             col_variable, 
                             p_title,
                             is_first){
  text_title = ""
  if(is_first){
    text_title = p_title
  }
  
  tmp_position = 'stack'
  tmp_y_label = 'Count'
  if(show_freq){
    tmp_position = 'fill'
    tmp_y_label = '%'
  }
  
  
  plt <- ggplot(data = p_data_plot, aes(x=id_x, y=count, fill=factor(col)))+
    geom_bar(stat='identity', position=tmp_position)+
    theme_minimal()+
    theme(legend.title = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.text.x=element_text(angle = 90, hjust = 0, vjust = 0.5),
          axis.title.x = element_blank())+ 
    #panel.grid.minor = element_blank())+
    scale_fill_manual(values=colors_fixed)+
    #xlab(x_variable)+
    ylab(tmp_y_label)+
    ggtitle(text_title)+
    scale_x_continuous(breaks=1:length(x_axis),
                       labels=factor(x_axis),
                       limits=c(0,length(x_axis)+1))
  
  
  return(plt)
}