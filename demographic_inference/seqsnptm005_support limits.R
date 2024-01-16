# Demographic parameter estimation was done for a grid of parameter combinations and then
# interpolated (separate Mathematica script). This script gets the support limits for the
# demographic parameter estimates by taking the quantiles of 500,000 draws from the interpolated
# likelihood surface.
#
#author: "Anja Westram"
#
######################################################################################################

rm(list=ls())


# Get interpolation results & transform back
grid = read.table("../Data/Interpolation demographic parameters_random values 23 Sept 2023.csv",
                  header=F, stringsAsFactors=F, sep=",")
names(grid) = c("N0", "r", "K", "M", "nll")

grid$N0 = exp(grid$N0) # Starting population size
grid$r = exp(grid$r) # Growth parameter
grid$K = exp(grid$K) # Carrying capacity
grid$M = exp(grid$M) - 0.5 # Number of migrants



# Get means and quantiles
mean(grid$N0)
quantile(grid$N0, probs = 0.025)
quantile(grid$N0, probs = 0.975)
mean(grid$r)
quantile(grid$r, probs = 0.025)
quantile(grid$r, probs = 0.975)
mean(grid$K)
quantile(grid$K, probs = 0.025)
quantile(grid$K, probs = 0.975)
mean(grid$M)
quantile(grid$M, probs = 0.025)
quantile(grid$M, probs = 0.975)
