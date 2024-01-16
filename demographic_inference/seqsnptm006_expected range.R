# Get expected range of allele frequency change (under drift, gene flow and model uncertainty)
# for each SNP and inversion by simulating 1,000 replicates under parameter combinations drawn
# from the inferred likelihood surface. Then compare to observed allele frequency changes.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())

# Set which loci to work with
type = "INV" # "NEU", "OUT", "INV" for control SNPs, spatial outlier SNPs, inversions 



# Get random draws from likelihood surface of interpolation of demographic parameter likelihoods
grid = read.table("../Data/Interpolation demographic parameters_random values 23 Sept 2023.csv", header=F,
                  stringsAsFactors=F, sep=",")
names(grid) = c("N0", "r", "K", "M", "nll")
set.seed(5) # 5
grid = grid[sample(1:length(grid$nll), 1000, replace=F), ]



# Get observed allele frequency data
if(type=="NEU"|type=="OUT"){
  obs = read.table(paste("../Data/SEQSNPTM001_ALLELE_FREQS_", type, "_subsampled.txt", sep=""),
                 header=T, stringsAsFactors=F)
  obs$pC0000 = NULL # Crab population data not needed
}

if(type=="INV"){
  obs = read.table("../Data/Inversions_Frequencies_for_Envelope_November2023.txt",
                     header=T, stringsAsFactors=F)
  names(obs)[1] = "cp"
  obs = obs[, c("pC1992", "pW0000", "pS2005", "pS2018", "pS2021", "cp")]
  
  # Keep only one arrangement for complex inversions
  # (should be the one with the highest frequency in Wave)
  obs = obs[duplicated(obs$cp)==F, ]
}



# Polarise to make sure we always look at the allele with a higher frequency
# in Wave than Crab
rev = obs$pC1992 > obs$pW0000
obs[rev, 1:5] = 1 - obs[rev, 1:5]



# Set up data frame for results; each row is a locus; each column will give the simulated allele
# frequency at the end of the experiment in a single simulation replicate
res = as.data.frame(matrix(NA, nrow = length(obs$pC1992), ncol=1000))



# 1000 replicate simulations under neutrality
for (repl in 1:1000){
  print(repl)

  
  # Get parameter combinations from likelihood surface
  N0_set = exp(grid[repl, "N0"])
  r_set = exp(grid[repl, "r"])
  K_set = exp(grid[repl, "K"])
  M_set = exp(grid[repl, "M"]) - 0.5
  gen_factor = 2
  
  
  
  # Get total number of generations 
  gens = round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3)

  
  
  # Make population size for each generation using logistic growth function
  growth = function(gen) {K_set / (1 + ((K_set-N0_set)/N0_set) * exp(-r_set*gen))}
  Ns = round(growth(1:gens), 0)

  
  
  # Run allele frequency change for each locus under gene flow and drift
  for (locus in 1:length(obs$pC1992)){

    pW = obs$pW0000[locus] # Wave allele frequency for current locus
    pS = obs$pC1992[locus] # Skerry starting frequency for current locus (will change)
    N = N0_set # Starting population size
    
    
    
    # Run gene flow and drift for inferred number of generations
    for (gen in 1:gens){

      
      # Expected frequency after gene flow
      pp = ((N-M_set)*pS + M_set*pW) / N

      
      
      # Frequency after one generation of drift
      N = Ns[gen]
      pS = rbinom(1, N, pp) / N

    }
    
    
    
    # Add final allele frequency on the skerry from this simulation run to results data frame
    if(type=="OUT"|type=="NEU"){
      res[locus, repl] = rbinom(1, 68, pS)/68 # Sampling of 34 individuals
    }
    
    if(type=="INV"){ # Consider sample size variation between inversions
      ns = read.table("../Data/Inversions_SampleSize_for_Envelope.txt", header=T,
                      stringsAsFactors=F, sep="\t") # Diploid sample sizes
      n = ns$n[ns$Inversion==obs$cp[locus] & ns$Population=="Skerry 2021"] * 2 # Haploid sample size
                                                                               # for current inversion
      res[locus, repl] = rbinom(1, n, pS)/n
    }
  }
}



# Quantiles and median of expected allele frequency at the end of the experiment (2021)
obs$ql = apply(res, 1, function(x) quantile(x, 0.025))
obs$qh = apply(res, 1, function(x) quantile(x, 0.975))
obs$med = apply(res, 1, function(x) median(x))



# Quantiles and median of expected allele frequency *change*
obs$afd = obs$pS2021-obs$pC1992
obs$afdl = obs$ql - obs$pC1992
obs$afdh = obs$qh - obs$pC1992
obs$afdm = obs$med - obs$pC1992



# Plot observed allele frequency change and expected range for each locus, order by obs. allele
# frequency change. Expected range in dark blue if outside expected range and below median,
# light blue if inside expected range and below median, light red if inside expected range
# and above median, dark red if outside expected range and above median.
# (Export 8x4)
obs = obs[order(obs$afd), ]
plot(obs$afd, xaxt="n", ylab="Allele frequency change on skerry", ylim=c(-.5,1), xlab="", cex=.5)
for(num in 1:length(obs$pC1992)){ 
  afd  = obs$afd[num]
  afdl = obs$afdl[num]
  afdh = obs$afdh[num]
  afdm = obs$afdm[num]
  if(afd>afdl & afd<afdh & afd<=afdm){lines(c(num,num), c(afdl,afdh), lwd=1, col=rgb(0.68,0.85,1))}
  if(afd>afdl & afd<afdh & afd>afdm){lines(c(num,num), c(afdl,afdh), lwd=1, col=rgb(1,0.7,0.7))}
  if(afd<=afdl){lines(c(num,num), c(afdl,afdh), lwd=1, col=rgb(0.2,0.2,1))}
  if(afd>=afdh){lines(c(num,num), c(afdl,afdh), lwd=1, col=rgb(0.8,0,0))}
}

points(obs$afd, cex=.5, col="black", pch=16) # Add observation on top again
lines(c(-1000,1000), c(0,0)) # Line at 0 change

if(type=="INV"){ # Label with inversion names
  axis(1, 1:length(obs$cp), obs$cp, las=2)
}



# Proportions above and below quantiles and median
length(obs$pS2021[obs$pS2021>obs$ql & obs$pS2021<obs$qh]) / length(obs$pS2021)
length(obs$pS2021[obs$pS2021>=obs$qh]) / length(obs$pS2021)
length(obs$pS2021[obs$pS2021<=obs$ql]) / length(obs$pS2021)
length(obs$pS2021[obs$pS2021>obs$med]) / length(obs$pS2021)



# Write out which loci are above expected range
obs$out = obs$afd>=obs$afdh
write.table(obs[,c("cp","out")], paste("../Data/SEQSNPTM006_OUTLIER_STATUS_", type, ".txt", sep=""),
            quote=F, row.names=rownames(obs), col.names=F)
