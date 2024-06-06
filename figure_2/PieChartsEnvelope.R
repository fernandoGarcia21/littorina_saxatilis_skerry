############################## Pie Charts Figure 2A ############################
# Creates two pie charts: control and spatial outliers
# The figure shows the percentage of collinear loci beyond 
# and within neutral expectation based on the demographic inference
#
#author: "Diego Garcia"
#date: "2023-05-11"
################################################################################

# Percentages of SNPs with AF change above median envelope
perc_neutral_above = 0.55
#perc_outlier_above = 0.67 #without downsampling of outliers to reduce clustering (full spatial outlier dataset)
perc_outlier_above = 0.62 #with downsampling of outliers to reduce clustering (reduced-LD spatial outlier dataset)

# Percentages of control SNPs with AF change beyond the expected range
perc_neutral_above_median_beyond_exp = 0.05
perc_neutral_below_median_beyond_exp = 0.06

# Percentages of spatial outlier SNPs with AF change beyond the expected range
#perc_outlier_above_median_beyond_exp = 0.09 #without downsampling of outliers to reduce clustering
#perc_outlier_below_median_beyond_exp = 0.01 #without downsampling of outliers to reduce clustering
perc_outlier_above_median_beyond_exp = 0.14 #without downsampling of outliers to reduce clustering
perc_outlier_below_median_beyond_exp = 0.04 #without downsampling of outliers to reduce clustering

#Colour scheme of the pie charts
pies_palette = c('Above median within range' = '#73BADA',
                 'Above median beyond range' = '#4674B7', 
                 'Below median within range' = '#F5F5DC',
                 'Below median beyond range' = '#F5CE79')

#Data frame with the information (percentages) of control loci
neutral_pie_df <- data.frame(perc=c(perc_neutral_above - perc_neutral_above_median_beyond_exp, 
                                    perc_neutral_above_median_beyond_exp,
                                    (1 - perc_neutral_above) - perc_neutral_below_median_beyond_exp,
                                    perc_neutral_below_median_beyond_exp),
                             serie=c('Above median within range',
                                     'Above median beyond range',
                                     'Below median within range',
                                     'Below median beyond range'))

#Data frame with the information (percentages) of spatial outlier loci
outlier_pie_df <- data.frame(perc=c(perc_outlier_above - perc_outlier_above_median_beyond_exp, 
                                    perc_outlier_above_median_beyond_exp,
                                    (1 - perc_outlier_above) - perc_outlier_below_median_beyond_exp,
                                    perc_outlier_below_median_beyond_exp),
                             serie=c('Above median within range',
                                     'Above median beyond range',
                                     'Below median within range',
                                     'Below median beyond range'))


###################################################################
# Returns a grid with the control and outlier pie charts 
# from the neutral envelope statistics
###################################################################
generate_envelope_pie_charts <- function(p_show_title){
  
  plot_pie_neutral <- plot_pie(neutral_pie_df, 'Control', pies_palette, FALSE, 'Above median beyond range', p_show_title)
  plot_pie_outlier <- plot_pie(outlier_pie_df, 'Spatial outliers', pies_palette, TRUE, 'Above median beyond range', p_show_title)
  
  legend_outlier = cowplot::get_legend(plot_pie_outlier)
  
  row_pies = grid.arrange(plot_pie_neutral, 
                          plot_pie_outlier+theme(legend.position = 'none'), 
                          nrow = 1, ncol = 2)
  
  grid_pies = grid.arrange(row_pies, 
                           legend_outlier,
                           nrow = 2, ncol = 1,
                           heights=c(0.9,0.1))
  return(grid_pies)
}


###################################################################
# Generates a pie chart according to the percentages of each series.
# p_right_side_serie is the name of the series that is centered at the most
# right side of the pie chart. 
###################################################################
plot_pie <- function(data_df, p_title, p_palette, is_second_pie, p_right_side_serie, p_show_title){
  
  #Determine the coordinate in radians for starting to plot the pie chart
  tmp_rad_start = estimate_rad_start_drawing(data_df, p_right_side_serie)
  aux_top_margin = -1
  if(p_show_title){
    aux_top_margin = 0
  }
  
  #Generate a ggplot object (pie chart)
  plt_pie <- ggplot(data_df, aes(x = "", y = perc, fill = factor(serie))) +
    geom_col(width = 1, size = 0.0, color = "white") +
    geom_text(aes(label = paste0(round(perc*100, digits = 1),'%')),
              position = position_stack(vjust = 0.5),
              show.legend = FALSE,
              color = 'white',
              size=4
    ) +
    coord_polar(theta = "y", start=tmp_rad_start)+
    scale_fill_manual(values = p_palette)+
    theme_void()+
    ggtitle(p_title)+
    theme(legend.title = element_blank(),
          plot.margin = margin(aux_top_margin,0,0,0, "cm"))
  
  if(is_second_pie){
    plt_pie <- plt_pie+theme(legend.position = 'bottom')
  }else{
    plt_pie <- plt_pie+theme(legend.position = 'none')
  }
  
  return(plt_pie)
  
}


###################################################################
# Estimate the coordinate in radians where the circle starts.
# p_right_side_serie is the name of the series that is centered at the most
# right side of the pie chart. It is used as a reference to estimate
# the length of that segment and then where to start the circle.
###################################################################
estimate_rad_start_drawing <- function(data_df, p_right_side_serie){
  #Estimate the length of the circunference in rads
  tmp_rad_circunf = (2*pi)
  
  perc_right = data_df[data_df$serie == p_right_side_serie, ]$perc
  #Estimate the length of the segment of outside neutral expectation, in rads
  tmp_rad_outneutral = tmp_rad_circunf*perc_right
  
  #The figure starts 
  if(perc_right > 0.70){
    tmp_start_drawing = (tmp_rad_circunf*0.75)+((tmp_rad_circunf-tmp_rad_outneutral)/2)
  }else{
    tmp_start_drawing = (tmp_rad_circunf*0.25)+(tmp_rad_outneutral/2)
  }
  
  
  return(tmp_start_drawing)
}


