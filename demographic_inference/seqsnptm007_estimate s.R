# Estimate selection coefficients for spatial outliers and inversions above the expected range.
# For that, use transition matrix similar to demographic inference, but fix demographic
# parameters to the maximum-likelihood ones.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())

library(Matrix)


# Set which loci to work with
type="inv" # "out" for spatial outliers, "inv" for inversions



# Set demographic parameters
N0_set = 54 # Starting population size
r_set = 0.12 # Growth parameter
K_set = 1329 # Carrying capacity
M_set = 3.24 # Number of migrants
gen_factor = 2 # Generations per year
gens = round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3) # Total number of generations
sampling_gens = c(round(gen_factor*13),
                  round(gen_factor*13)+round(gen_factor*13),
                  round(gen_factor*13)+round(gen_factor*13)+round(gen_factor*3)) # Samplin generations
s_list = seq(0,2,0.02) # Values of selection coefficient to test



# Get observed allele frequencies (Crab, Wave and 3 years for Skerry) and keep only loci
# above expected range
if(type=="out"){
  
  # Get observed allele frequencies and keep only loci with evidence for selection
  obs = read.table("../Data/SEQSNPTM001_ALLELE_FREQS_OUT_subsampled.txt", 
                   header=T, stringsAsFactors=F)
  outlier_status = read.table("../Data/SEQSNPTM006_OUTLIER_STATUS_OUT.txt", header=F,
                              stringsAsFactors = F) # Whether each locus is above expected range
  names(outlier_status) = c("num", "cp", "outlier_status")
  
  obs = obs[obs$cp %in% outlier_status$cp[outlier_status$outlier_status==T], ]
}

if(type=="inv"){
  
  # Get observed allele frequencies and rearrange, keep only loci with evidence for selection
  obs = read.table("../Data/Inversions_Frequencies_for_Envelope_November2023.txt",
                   header=T, stringsAsFactors=F, sep="\t")
  names(obs)[1] = "cp"
  obs = obs[, c("pC1992", "pW0000", "pS2005", "pS2018", "pS2021", "cp")]
  obs = obs[duplicated(obs$cp)==F, ] # Keep only arrangement more common in Wave, which is the first one
  
  outlier_status = read.table("../Data/SEQSNPTM006_OUTLIER_STATUS_INV.txt", header=F,
                              stringsAsFactors = F) # Whether each locus is above expected range
  names(outlier_status) = c("num", "cp", "outlier_status")
  
  obs = obs[obs$cp %in% outlier_status$cp[outlier_status$outlier_status==T], ]
}



# Polarise to make sure we always look at the allele with a higher frequency
# in Wave than Crab
rev = obs$pC1992 > obs$pW0000
obs[rev, 1:5] = 1 - obs[rev, 1:5]



# Make columns that will contain nlls for different selection coefficients
obs[,paste("s_set", s_list, sep="")] = NA



#############################################################################################################
# FUNCTIONS & POP GROWTH ####################################################################################
#############################################################################################################

# Function to get probability of sampling allele A after gene flow & selection
gf = function(p){
  pp = ((oldN-M_set)*p + M_set*pW) / oldN # pW is the Wave allele freq for the current locus
  pp = pp*(1+s) / (pp*s+1) # Selection
  return(pp)
}



# Make population size for each generation
growth = function(gen) {K_set / (1 + ((K_set-N0_set)/N0_set) * exp(-r_set*gen))}
Ns = round(growth(1:gens), 0)
plot(1:gens, Ns, ylim=c(0,K_set+1000))



# Get time when the carrying capacity is reached, as from then on we do not need to re-calculate the
# matrix in every generation again
carrying = which(Ns==K_set)[3] # Generation from which TM is stable
carrying[is.na(carrying)] = 99999999 # Set to arbitrarily high value if carrying capacity is never reached



#############################################################################################################
# RUN FOR EACH LOCUS ###########################################################################################
#############################################################################################################

plot(0,0, xlim=c(0,2), ylim=c(0,1000), col="white") # Plot that will contain nll for each SNP & s

for (locus in 1:length(obs$pC)){
  
  for (s in s_list){ # Loop through selection coefficients
    
    print(paste(locus, s))
    
    pW = obs$pW0000[locus] # Allele freq in Wave for this locus
    
    
    
    # Obtain transition matrices for the three time intervals (1992-2005, 2005-2018, 2018-2021)
    M = Matrix(diag(N0_set+1), sparse=T) # Identity matrix to multiply with in the first generation
    for (gen in 1:gens){ 
      
      
      
      # Get previous and current population size
      if(gen==1){oldN=N0_set}else{oldN=Ns[gen-1]}
      newN = Ns[gen]
      
      
      
      # Get TM describing change by drift and gene flow in 1 generation,
      # i.e. probability of going from i to j copies.
      # i is the old number of copies, ranging between 0 and oldN.
      # j is the new number of copies, ranging between 0 and newN.
      # Value is determined by binomial sampling of newN copies from the old frequency.
      if(gen<carrying){ # Generate new matrix to multiply with only as long as pop size changes)
        #print(paste(newN, K_set))
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
    if(type=="out"){
      s1992 = dbinom(obs$pC1992[locus]*36, 36, seq(0,1,length=dim(M1992_2005)[1]))
      s2005 = dbinom(obs$pS2005[locus]*62, 62, seq(0,1,length=dim(M1992_2005)[2]))
      s2018 = dbinom(obs$pS2018[locus]*90, 90, seq(0,1,length=dim(M2005_2018)[2]))
      s2021 = dbinom(obs$pS2021[locus]*68, 68, seq(0,1,length=dim(M2018_2021)[2]))
    }
    
    if(type=="inv"){ # (Use relevant sample size for each inversion)
      
      ns = read.table("../Data/Inversions_SampleSize_for_Envelope.txt", header=T,
                      stringsAsFactors=F, sep="\t") # Diploid sample sizes
      n1992 = ns$n[ns$Inversion==obs$cp[locus] & ns$Population=="Crab 1992"] * 2 # Haploid sample size for this inversion 
      n2005 = ns$n[ns$Inversion==obs$cp[locus] & ns$Population=="Skerry 2005"] * 2
      n2018 = ns$n[ns$Inversion==obs$cp[locus] & ns$Population=="Skerry 2018"] * 2
      n2021 = ns$n[ns$Inversion==obs$cp[locus] & ns$Population=="Skerry 2021"] * 2
      
      s1992 = dbinom(obs$pC1992[locus]*n1992, n1992, seq(0,1,length=dim(M1992_2005)[1]))
      s2005 = dbinom(obs$pS2005[locus]*n2005, n2005, seq(0,1,length=dim(M1992_2005)[2]))
      s2018 = dbinom(obs$pS2018[locus]*n2018, n2018, seq(0,1,length=dim(M2005_2018)[2]))
      s2021 = dbinom(obs$pS2021[locus]*n2021, n2021, seq(0,1,length=dim(M2018_2021)[2]))
    }
    
    
    
    # For each time interval, get the probability of going from i to j copies under gene flow and drift
    # and then sampling the observed count
    ms2005 = t(t(M1992_2005) * s2005) # Simple product!
    ms2018 = t(t(M2005_2018) * s2018)
    ms2021 = t(t(M2018_2021) * s2021)
    
    
    
    # Get nll - based on all intervals
    # First part is uniform prior for true allele frequencies in 1992
    obs[,paste("s_set", s, sep="")][locus] = -log(sum((rep(1/(N0_set+1), (N0_set+1)) * s1992) %*% ms2005 %*% ms2018 %*% ms2021)) 
  }
  
  
  
  lines(s_list, obs[locus,paste("s_set", s_list, sep="")],
        col=rgb(0.5,0.5,0.5,0.2)) # Plot nlls for different s for this SNP
  
}



# Write out
write.table(obs, paste("../Data/SEQSNPTM007_ESTIMATE_S_", type, ".txt", sep=""),
            quote=F, col.names=T, row.names=T)