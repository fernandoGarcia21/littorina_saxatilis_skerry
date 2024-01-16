#
#
#  run the Stan model for selection using the Skerry data
#
#  run is for phenotype 'phen' which is selected in the <read_data.R> script
#
#  @author: Roger Butlin
#




# recommended rstan start-up options
library("rstan") 
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# add loo and bayesplot
library("loo")
library("bayesplot")

setwd(".")  # directory containing the Stan file and data files

nchains <- 4

gen_var <- Vg_site[phen] # choose one of these to set the Vg assumption
gen_var <- Vg_crab[phen]

# set up input data for the Stan call

stan_data <- list(g = 2, # set the generations per year assumption here
                  Ng = 29,
                  Nd = length(donor$snailID),
                  Ns96 = length(skerry96$snailID),
                  Ns02 = length(skerry02$snailID),
                  Ns05 = length(skerry05$snailID),
                  Ns18 = length(skerry18$snailID),
                  Ns21 = length(skerry21$snailID),
                  Nw = length(wave$snailID),
                  don = donor[,phen],  
                  S96 =  skerry96[,phen], 
                  S02 =  skerry02[,phen], 
                  S05 =  skerry05[,phen], 
                  S18 =  skerry18[,phen], 
                  S21 =  skerry21[,phen], 
                  wav = wave[,phen],   
                  p_est =  as.numeric(p[phen]),
                  sdp =  as.numeric(sd_p[phen]),
                  Vg_est = as.numeric(gen_var), 
                  sd_Vg = as.numeric(sd_Vgs[phen]))

# set starting values for the chains

init_fun <- function(chainID=1) list(p=as.numeric(p[phen]),
                                     Vg=as.numeric(gen_var),
                                     Vs=20*as.numeric(gen_var),
                                     crab=mean(donor[,phen]),
                                     wave_opt=mean(wave[,phen]),
                                     Ve=4*as.numeric(gen_var))
init_ll <- lapply(1:nchains, function(id) init_fun(chainID = id))

###### run Stan fit (NB flat prior on Vs/Vp from 0 to 50)
fit_Vgs <- stan(
  file = "skerry_g_p.stan",  # Stan program - this version sets an upper bound on Vs/Vp
  data = stan_data,      # named list of data
  chains = nchains,             # number of Markov chains
  warmup = 1000,          # number of warmup iterations per chain
  iter = 3000,            # total number of iterations per chain
  cores = 4,              # number of cores (could use one per chain)
  refresh = 500,          # progress shown
  control = list(adapt_delta = 0.95)
)

# print summary of fit
fit_Vgs_summary <- summary(fit_Vgs,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975)) 
fit_Vgs_summary$summary[c("crab", "wave_opt", "p", "Vg","Vs","sd_phen","Vs_Vp","w_crab"),]

# get values to plot

fit_mean <- summary(fit_Vgs, pars = c("fit_phen"))$summary[,'mean']
fit_uci <- summary(fit_Vgs, pars = c("fit_phen"))$summary[,'97.5%']
fit_lci <- summary(fit_Vgs, pars = c("fit_phen"))$summary[,'2.5%']

fit_year <- 1991 + seq(1,stan_data[[2]]+1,1/stan_data[[1]])

plotdata <- data.frame(year,means,upper,lower)
fitdata <- data.frame(year=fit_year,means=fit_mean,fit_lci,fit_uci)

# set colours
cols = c('black', rep('red', 5), 'black')

# plot
ggplot(data = fitdata, aes(x=year,y=means), show.legend =
         FALSE) +
  # line inherits the x and y and data; just specify colour
  geom_line(colour = 'red') +
  # ribbon inherits x and y and data; just specify new aes
  geom_ribbon(aes(ymax = fit_uci, ymin = fit_lci), fill = 'pink', alpha = 0.5) +
  # new data to use here....
  # because the names are the same, need only specify the data and the cols
  geom_point(data = plotdata, col = cols) +
  # new data here, and new aes()
  geom_errorbar(data = plotdata, aes(ymin = lower, ymax =
                                       upper), col = cols)+
  # adjust the labels
  labs(x = "Year", y = phen)+
  # old skool and bigger fonts.
  theme_classic(base_size = 15)


  
