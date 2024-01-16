# Calculate the total negative log-likelihood for each parameter combination used in the demographic
# inference by summing over SNPs.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())


# Get all parameter combinations
grid = read.table("../Data/SEQSNPTM002_GRID.txt", header=F, stringsAsFactors=F)
names(grid) = c("N0_set", "r_set", "K_set", "M_set", "gen_factor")
grid$nll = NA



# Get nll for each parameter combination (transition matrix results) and SNP, and sum over
# SNPs to get a single value per parameter combination
for (lin in 1:length(grid$N0_set)){
  print(lin)
  
  N0_set = grid$N0_set[lin]
  r_set = grid$r_set[lin]
  K_set = grid$K_set[lin]
  M_set = grid$M_set[lin]
  gen_factor = grid$gen_factor[lin]
  
  fil = paste("../Data/SEQSNPTM003_DEMOGRAPHIC INFERENCE/SEQSNPTM003_nll", N0_set, r_set, K_set, round(M_set,2), round(gen_factor, 2), ".txt", sep="_")
  
  if (file.exists(fil)){
    
    dat = read.table(paste("../Data/SEQSNPTM003_DEMOGRAPHIC INFERENCE/SEQSNPTM003_nll", N0_set, r_set, K_set, round(M_set,2), round(gen_factor,2), ".txt", sep="_"))
    
    grid[lin, "nll"] = sum(dat$nll)
  }
}



# Check for completeness
length(grid$nll[is.na(grid$nll)])



# Check best model
grid[grid$nll==min(grid$nll, na.rm=T) & is.na(grid$nll)==F, ]



# Write out
write.table(grid, "../Data/SEQSNPTM004_RESULTS DEMOGRAPHIC INFERENCE.txt", col.names=T, row.names=F, quote=F)
