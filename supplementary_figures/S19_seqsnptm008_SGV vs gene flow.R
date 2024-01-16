# Simulate allele frequency change at each spatial outlier and inversion outside the expected
# range (i.e. with evidence for selection on the skerry) under the previously inferred demogra-
# phic parameters and s. Simulate allelic copies from standing genetic variation (SGV) and gene
# flow (GF) separately to estimate the contribution of the two to adaptation.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())


# Type of locus to analyse
type = "inv" # "out" for spatial outliers, "inv" for inversions



# Get best estimate of s for each locus from transition matrix inference
if(type == "out"){
  obs = read.table("../Data/SEQSNPTM007_ESTIMATE_S_out.txt", header=T,
                   stringsAsFactors = F)
}
if(type == "inv"){
  obs = read.table("../Data/SEQSNPTM007_ESTIMATE_S_inv.txt", header=T,
                   stringsAsFactors = F)
}
obs$best_s = apply(obs[, 7:107], 1, function(x) names(which(x==min(x))))
obs$best_s = as.numeric(gsub("s_set", "", obs$best_s))



# Get best estimate of s for each locus based on interpolation of TM results,
# and compare to TM best estimate 
if (type == "out"){
  interp = read.table("../Data/Interpolation s sel SNPs 10 Dec.csv",
                      sep=",", header=F, stringsAsFactors=F)
}
if (type == "inv"){
  interp = read.table("../Data/Interpolation s inversions 10 Dec.csv",
                      sep=",", header=F, stringsAsFactors=F)
}
names(interp) = c("s_min", "s", "s_max")
interp$cp = obs$cp

plot(obs$best_s, interp$s)



# Get observed allele frequency data and best s from interpolation into same data frame
obs = merge(obs, interp, by="cp")



# Inferred demographic parameters
N0_set = 54
r_set = 0.12
K_set = 1329
M_set = 3.24
gen_factor = 2
gens = round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3)
sampling_gens = c(round(gen_factor*13), # Generations for which we have allele frequencies
                  round(gen_factor*13)+round(gen_factor*13),
                  round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3))



#############################################################################################################
# Simulate contributions of SGV and GF ######################################################################
#############################################################################################################

# Function to get probability of sampling allelic copies after selection
sel = function(p0, p1, p2){ # p0 = Crab allele; p1 = Wave allele from SGV; p2 = Wave allele from gene flow
  pp0 = p0*1     / (p0*1 + p1*(1+s) + p2*(1+s)) # Selection
  pp1 = p1*(1+s) / (p0*1 + p1*(1+s) + p2*(1+s)) # Selection
  pp2 = p2*(1+s) / (p0*1 + p1*(1+s) + p2*(1+s)) # Selection
  return(c(pp0, pp1, pp2))
}



# Make population size for each generation
growth = function(gen) {K_set / (1 + ((K_set-N0_set)/N0_set) * exp(-r_set*gen))}
Ns = round(growth(1:gens), 0)
plot(1:gens, Ns, ylim=c(0,K_set+1000))



# Prepare plot
pdf(paste("../Data/SEQSNPTM008_SGV_VS_GF_", type, ".pdf", sep=""), 10, 16)
par(mfrow=c(7,3), mar=c(5,4,2,2), oma=c(3,3,3,3))



# Run 1000 simulations under gene flow, selection and drift for each locus, modelling allelic
# copies from gene flow vs. standing genetic variation separately
for (snp in 1:length(obs$cp)){
  
  
  # Allele freq in Wave for this SNP
  pW = obs$pW0000[snp]
  
  
  # s estimated for this SNP
  s = obs$s[snp]
  
  
  
  # Matrices that will contain simulated allele frequencies for this SNP for each generation and replicate
  res0 = matrix(ncol=gens, nrow=1000) # Allele frequency of Crab allele
  res1 = matrix(ncol=gens, nrow=1000) # Allele frequency of Wave allele from SGV
  res2 = matrix(ncol=gens, nrow=1000) # Allele frequency of Wave allele from gene flow
  
  
  
  # Run 1000 replicates
  for(repl in 1:1000){
    
    
    # Starting allele frequencies
    p0 = 1-obs$pC1992[snp] # Allele frequency of Crab allele
    p1 = obs$pC1992[snp] # Allele frequency of Wave allele from SGV
    p2 = 0 # Allele frequency of Wave allele from gene flow
    
    
    
    # Run for all generations of experiment
    for (gen in 1:gens){ 
      
      
      
      # Get previous and current population size
      if(gen==1){oldN=N0_set}else{oldN=Ns[gen-1]}
      newN = Ns[gen]
      
      
      
      # Expected frequencies after gene flow...
      p0 = ((oldN - M_set) * p0 + (M_set * (1-pW))) / oldN
      p1 = ((oldN - M_set) * p1) / oldN
      p2 = ((oldN - M_set) * p2 + M_set * pW) / oldN
      
      
      
      # ... and selection
      pp0 = p0*1     / (p0*1 + p1*(1+s) + p2*(1+s))
      pp1 = p1*(1+s) / (p0*1 + p1*(1+s) + p2*(1+s))
      pp2 = p2*(1+s) / (p0*1 + p1*(1+s) + p2*(1+s))
      
      
      
      # New frequencies after drift
      alleles = sample(c(0,1,2), newN, T, prob=c(pp0,pp1,pp2))
      p0 = length(alleles[alleles==0]) / length(alleles)
      p1 = length(alleles[alleles==1]) / length(alleles)
      p2 = length(alleles[alleles==2]) / length(alleles)
      
      
      
      # Add to output     
      res0[repl, gen] = p0
      res1[repl, gen] = p1
      res2[repl, gen] = p2
      
    }
    
  }
  
  
  
  # Get average allele frequencies across simulations
  sim0 = colMeans(res0)
  sim1 = colMeans(res1)
  sim2 = colMeans(res2)
  
  
  
  # Plot total simulated Wave allele frequency trajectory, averaged over replicates
  plot(0:gens, c(obs$pC1992[snp], sim1 + sim2), col="black", ylim=c(0,1), type="l", lwd=2,
       main=obs$cp[snp], xlab="Generation", ylab="Allele frequency", yaxt="n") # Wave allele frequency from SGV
  axis(2, seq(0,1,0.2), seq(0,1,0.2), las=2)
  
  
  
  # Add simulated contribution of SGV vs gene flow, from simulations and from back-
  # of-envelope approximation, except for loci where observations and simulations differ
  # a lot
  if(obs$cp[snp] %in% c("Contig91839_37927", "Contig55826_16617")==F){
    lines(1:gens, sim2, col=rgb(115,186,218,m=255), lwd=2, lty=2) # Wave allele frequency from gene flow
    lines(1:gens, sim1, col=rgb(224,145,87,m=255), lwd=3, lty=2) # Total Wave allele frequency from SGV
    
    
    
    # Back-of-envelope expectation - ratio (See Supplementary Text for derivation)
    sgv = obs$pC1992[snp]
    gf = (M_set/N0_set) * obs$pW0000[snp] / (r_set + obs$s[snp])
    
    
    
    # Back-of-envelope expectation - expected frequency
    sgv_exp = sgv / (sgv+gf) * obs$pS2021[snp] # (sim1[58] + sim2[58])
    gf_exp = gf / (sgv+gf) * obs$pS2021[snp] # (sim1[58] + sim2[58])
    
    
    
    # Add to plot
    points(sampling_gens[3], sgv_exp, pch=16, cex=1.5, col=rgb(224,145,87,170,m=255))
    points(sampling_gens[3], gf_exp, pch=16, cex=1.5, col=rgb(115,186,218,170,m=255))
    
  }
  
  
  
  # Add observed allele frequencies
  points(0, obs$pC1992[snp], pch=16, cex=1.5)
  points(sampling_gens[1], obs$pS2005[snp], pch=16, cex=1.5)
  points(sampling_gens[2], obs$pS2018[snp], pch=16, cex=1.5)
  points(sampling_gens[3], obs$pS2021[snp], pch=16, cex=1.5)
  
  
  
  # Include estimated s into the plot
  mtext(paste("s=", round(obs$s[snp],2), sep=""), 1, line = -2)
  
}
dev.off()
