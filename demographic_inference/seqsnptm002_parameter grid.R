# Prepare a grid of all parameter combinations to test in the demographic inference.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())

N0s = c(20, 40, 80, 160, 320) # Starting haploid population size
rs = c(0.025, 0.05, 0.1, 0.2, 0.4) # Growth parameter
Ks = c(250, 500, 1000, 2000, 4000) # Carrying capacity
Ms = c(0, 1, 2, 4, 8) # Numbers of haploid migrants
gen_factors = c(4/3, 5/3, 2) # Number of generations per year

grid = expand.grid(N0s, rs, Ks, Ms, gen_factors) # All combinations

names(grid) = c("N0", "r", "K", "M", "gen_factor") # Add header

write.table(grid, "../Data/SEQSNPTM002_GRID.txt", col.names=F, row.names=F, quote=F) # Write out
