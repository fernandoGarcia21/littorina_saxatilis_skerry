# Main demographic inference script to get maximum likelihood for different combinations
# of demographic parameters, based on the allele frequency changes at control loci in the
# skerry over time, calculating allele frequency changes using transition matrices. This script
# runs the analysis for a single parameter combination and outputs the negative log-likelihood
# for each SNP.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())

library(Matrix)



# Set parameters manually:
N0_set = 80
r_set = 2
K_set = 100
M_set = 3
gen_factor = 4/3



# Alternatively, get parameters from wrapper (script runs a single combination)
# args = commandArgs(trailingOnly=TRUE)
# N0_set = as.numeric(args[1])     # Starting haploid Ne
# r_set = as.numeric(args[2])      # Growth parameter
# K_set = as.numeric(args[3])      # Carrying capacity
# M_set = as.numeric(args[4])      # Number (not rate) of haploid migrants per generation
# gen_factor = as.numeric(args[5]) # Number of generations 1 year corresponds to



# Calculated parameters
gens = round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3) # Total number of generations
sampling_gens = c(round(gen_factor*13), # Generations for which we have data
                  round(gen_factor*13)+round(gen_factor*13),
                  round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3))



# Get observed allele frequencies at control SNPs (Crab, Wave and 3 years for Skerry)
obs = read.table("../Data/SEQSNPTM001_ALLELE_FREQS_NEU_subsampled.txt",  header=T, stringsAsFactors=F)



# Use only 1 SNP per contig to reduce probability of LD
contigs = sapply(strsplit(obs$cp, "_", fixed=T), function(x) x[[1]])
obs = obs[duplicated(contigs)==F, ]



# Add column for negative log-likelihood of observed Skerry frequency under the current
# parameter combination
obs$nll = NA 



#############################################################################################################
# FUNCTIONS & POPULATION GROWTH #############################################################################
#############################################################################################################

# Function to get probability of sampling allele A after one generation of gene flow from Wave,
# given current frequency (p)
gf = function(p){
  pp = ((oldN-M_set)*p + M_set*pW) / oldN # pW is the Wave allele freq for the current locus
  return(pp)
}



# Make population size for each generation
growth = function(gen) {K_set / (1 + ((K_set-N0_set)/N0_set) * exp(-r_set*gen))} # Growth function
Ns = round(growth(1:gens), 0) # Population size for each generation
#plot(1:gens, Ns, ylim=c(0,K_set+1000))



# Get time when the carrying capacity is reached, as from then on we do not need to re-calculate the
# matrix in every generation again
carrying = which(Ns==K_set)[3] # Generation from which TM is stable
carrying[is.na(carrying)] = 99999999 # Set to arbitrarily high value if carrying capacity is never reached



#############################################################################################################
# RUN CALCULATION OF NEGATIVE LOG-LIKELIHOOD FOR EACH SNP ###################################################
#############################################################################################################

for (snp in 1:length(obs$pC)){
  
  
  
  # Allele frequency in Wave for this SNP
  pW = obs$pW0000[snp]
  
  
  
  # Obtain transition matrices for the three time intervals (1992-2005, 2005-2018, 2018-2021)
  M = Matrix(diag(N0_set+1), sparse=T) # Identity matrix to multiply with in the first generation
  for (gen in 1:gens){ 
    
    
    
    # Get previous and current population size
    if(gen==1){oldN=N0_set}else{oldN=Ns[gen-1]}
    newN = Ns[gen]
    
    
    
    # Get TM describing change by drift and gene flow in 1 generation,
    # i.e. probability of going from i (rows) to j (columns) copies.
    # i is the old number of copies, ranging between 0 and oldN.
    # j is the new number of copies, ranging between 0 and newN.
    # Value is determined by binomial sampling of newN copies from the old frequency.
    if(gen<carrying){ # Generate new matrix to multiply with only as long as pop size changes)
      M1 = t(sapply(0:oldN,
                    function(x) dbinom(0:newN, size=newN, prob=gf((x)/(oldN)))))
      
      # Remove very small values and use sparse matrix to speed up
      M1[M1<10^-9] = 0
      M1 = sweep(M1, 1, rowSums(M1), `/`)
      M1 = Matrix(M1, sparse=T)
    }
    
    
    
    M = M %*% M1 # Multiplication with matrix from previous generation - the matrix product
                 # across all generations gives the prob of going from i copies in the
                 # first to j copies in the last generation
    
    
    
    # If sampling happened in this generation, keep matrix and start new one if necessary
    if (gen == sampling_gens[1]){M1992_2005 = M; M=Matrix(diag(newN+1), sparse=T)}
    if (gen == sampling_gens[2]){M2005_2018 = M; M=Matrix(diag(newN+1), sparse=T)}
    if (gen == sampling_gens[3]){M2018_2021 = M}
  }
  
  
  
  # Vectors that give the probabilities of sampling the observed frequencies
  # at the 4 time points, under all possible true frequencies
  # (Numbers are haploid sample sizes)
  s1992 = dbinom(obs$pC1992[snp]*36, 36, seq(0,1,length=dim(M1992_2005)[1]))
  s2005 = dbinom(obs$pS2005[snp]*62, 62, seq(0,1,length=dim(M1992_2005)[2]))
  s2018 = dbinom(obs$pS2018[snp]*90, 90, seq(0,1,length=dim(M2005_2018)[2]))
  s2021 = dbinom(obs$pS2021[snp]*68, 68, seq(0,1,length=dim(M2018_2021)[2]))
  
  
  
  # For each time interval, get the probability of going from i to j copies under gene flow and drift
  # and then sampling the observed count
  ms2005 = t(t(M1992_2005) * s2005) # Simple product!
  ms2018 = t(t(M2005_2018) * s2018)
  ms2021 = t(t(M2018_2021) * s2021)
  
  
  
  # Get negative log-likelihood - based on all intervals
  # First part is uniform prior for true allele frequencies in 1992
  obs$nll[snp] = -log(sum((rep(1/(N0_set+1), (N0_set+1)) * s1992) %*% ms2005 %*% ms2018 %*% ms2021)) 

  
}



# Write out
write.table(obs, paste("../Data/SEQSNPTM003_DEMOGRAPHIC INFERENCE/SEQSNPTM003_nll", N0_set, r_set, K_set, round(M_set,2),
                       round(gen_factor,2), ".txt", sep="_"), quote=F)
