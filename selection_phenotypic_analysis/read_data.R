#
#
#
# set up data for the Stan fit to estimate selection on phenotypes on the Skerry
#
# @author: Roger Butlin
#
#


rm(list=ls())

options(stringsAsFactors = F)

library(plotly)

# data input and prep

setwd(".") # select data folder

donor1 <- read.csv("../Data/shape_parameters/CrabDonnor1992Parameters.csv")
donor2 <- read.csv("../Data/shape_parameters/DonnerPop1992_freezerParameters.csv")
donor3 <- read.csv("../Data/shape_parameters/DonnerPop2021Parameters.csv")
donor4 <- read.csv("../Data/shape_parameters/ExpSkerry2018_RamshC_Parameters.csv")

donor <- rbind(donor1,donor2,donor3,donor4)
donor$yr <- c(rep("92",21),rep("92f",25),rep("21",50),rep("18",33))

# check consistency among donor samples
pca <- princomp(donor[2:7])
plot(pca)
plot(pca$scores[donor$yr=="92",1],pca$scores[donor$yr=="92",2],xlim=c(-0.5,0.4),ylim=c(-0.3,0.4))
points(pca$scores[donor$yr=="92f",1],pca$scores[donor$yr=="92f",2],xlim=c(-0.5,0.4),ylim=c(-0.3,0.4),col="red")
points(pca$scores[donor$yr=="21",1],pca$scores[donor$yr=="21",2],xlim=c(-0.5,0.4),ylim=c(-0.3,0.4),col="blue")
points(pca$scores[donor$yr=="18",1],pca$scores[donor$yr=="18",2],xlim=c(-0.5,0.4),ylim=c(-0.3,0.4),col="green") 

skerry96 <- read.csv("../Data/shape_parameters/ExpSkerry1996Parameters.csv")
skerry02 <- read.csv("../Data/shape_parameters/ExpSkerry2002Parameters.csv")
skerry05 <- read.csv("../Data/shape_parameters/ExpSkerry2005Parameters.csv")
skerry21 <- read.csv("../Data/shape_parameters/ExpSkerry2021Parameters.csv")
skerry18 <- read.csv("../Data/shape_parameters/ExpSkerry2018_Rskerry_Parameters.csv")

wave18 <- read.csv("../Data/shape_parameters/ExpSkerry2018_SkareW_Parameters.csv")
wave21 <- read.csv("../Data/shape_parameters/ExpSkerry2021_SkareW_Parameters.csv")
wave <- rbind(wave18,wave21)

# rescale gw and gh to match parameters estimated by Koch et al (2022)
donor[,c("gw","gh")] <- exp(donor[,c("gw","gh")])
skerry96[,c("gw","gh")] <- exp(skerry96[,c("gw","gh")])
skerry02[,c("gw","gh")] <- exp(skerry02[,c("gw","gh")])
skerry05[,c("gw","gh")] <- exp(skerry05[,c("gw","gh")])
skerry18[,c("gw","gh")] <- exp(skerry18[,c("gw","gh")])
skerry21[,c("gw","gh")] <- exp(skerry21[,c("gw","gh")])
wave[,c("gw","gh")] <- exp(wave[,c("gw","gh")])

# if using log scale (for shellLength, gh, r0), convert to logs (i.e. reverse Koch et al. scaling for gh)
donor[,c("shellLength","gh","r0")] <- log(donor[,c("shellLength","gh","r0")])
skerry96[,c("shellLength","gh","r0")] <- log(skerry96[,c("shellLength","gh","r0")])
skerry02[,c("shellLength","gh","r0")] <- log(skerry02[,c("shellLength","gh","r0")])
skerry05[,c("shellLength","gh","r0")] <- log(skerry05[,c("shellLength","gh","r0")])
skerry18[,c("shellLength","gh","r0")] <- log(skerry18[,c("shellLength","gh","r0")])
skerry21[,c("shellLength","gh","r0")] <- log(skerry21[,c("shellLength","gh","r0")])
wave[,c("shellLength","gh","r0")] <- log(wave[,c("shellLength","gh","r0")])



# select phenotype for analysis
phen <- "c"


# summary plot of observations
year <- c(1992,1996,2002,2005,2018,2021,2023)

means <- c(mean(donor[,phen]),mean(skerry96[,phen]),mean(skerry02[,phen]),mean(skerry05[,phen]),mean(skerry18[,phen]),mean(skerry21[,phen]),mean(wave[,phen]))
vars <- c(var(donor[,phen]),var(skerry96[,phen]),var(skerry02[,phen]),var(skerry05[,phen]),var(skerry18[,phen]),var(skerry21[,phen]),var(wave[,phen]))
n <- c(length(donor[,phen]),length(skerry96[,phen]),length(skerry02[,phen]),length(skerry05[,phen]),length(skerry18[,phen]),length(skerry21[,phen]),length(wave[,phen]))
upper <- means+qt(0.975,n-1)*sqrt(vars/n)
lower <- means-qt(0.975,n-1)*sqrt(vars/n)

plotdata <- data.frame(year,means,upper, lower)

ggplot(plotdata, aes(year, means)) +        # ggplot2 plot with confidence intervals
  geom_point(colour=c("red",rep("black",5),"red")) +
  geom_errorbar(aes(ymin = lower, ymax = upper),colour=c("red",rep("black",5),"red"))


# get Koch et al. (2022) estimates of plasticity and genetic variance

setwd(".")  # data folder
CZA <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZA.txt",header = T)
CZB <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZB.txt",header = T)
CZD <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZD.txt",header = T)

# obtain mean estimates over transects and sd among transects for each phenotype
p <- (CZA[8,]+CZB[8,]+CZD[8,])/3 - (CZA[9,]+CZB[9,]+CZD[9,])/3  # plastic effect, unscaled
sd_p <- sqrt(((CZA[8,]-CZA[9,]-p)^2 + (CZB[8,]-CZB[9,]-p)^2 + (CZD[8,]-CZD[9,]-p)^2)/2) # among-site sd
Vg_site <- (CZA[5,]*CZA[13,]^2+CZB[5,]*CZB[13,]^2+CZD[5,]*CZD[13,]^2)/3  # need to multiply by s^2 to get back to original scale
sd_Vgs <- sqrt(((CZA[5,]*CZA[13,]^2-Vg_site)^2+(CZB[5,]*CZB[13,]^2-Vg_site)^2+(CZD[5,]*CZD[13,]^2-Vg_site)^2)/2)
Vg_crab <- (CZA[4,]*CZA[13,]^2+CZB[4,]*CZB[13,]^2+CZD[4,]*CZD[13,]^2)/3  # need to multiply by s^2 to get back to original scale
sd_Vgc <- sqrt(((CZA[4,]*CZA[13,]^2-Vg_crab)^2+(CZB[4,]*CZB[13,]^2-Vg_crab)^2+(CZD[4,]*CZD[13,]^2-Vg_crab)^2)/2)

params <- rbind(p, sd_p, Vg_site, sd_Vgs, Vg_crab, sd_Vgc) 

# or Koch et al. estimates on log scale (for traits analysed on this scale)

setwd(".")
CZA <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZA.log.txt",header = T)
CZB <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZB.log.txt",header = T)
CZD <- read.table("../Data/koch_2022_estimates/Crab.low.Wave.CZD.log.txt",header = T)

# obtain mean estimates over transects and sd among transects for each phenotype (log scale)
p <- (CZA[1,]+CZB[1,]+CZD[1,])/3 - (CZA[6,]+CZB[6,]+CZD[6,])/3  # plastic effect, unscaled
sd_p <- sqrt(((CZA[1,]-CZA[6,]-p)^2 + (CZB[1,]-CZB[6,]-p)^2 + (CZD[1,]-CZD[6,]-p)^2)/2) # among-site sd
Vg_site <- (CZA[5,]*CZA[10,]^2+CZB[5,]*CZB[10,]^2+CZD[5,]*CZD[10,]^2)/3  # need to multiply by s^2 to get back to original scale
sd_Vgs <- sqrt(((CZA[5,]*CZA[10,]^2-Vg_site)^2+(CZB[5,]*CZB[10,]^2-Vg_site)^2+(CZD[5,]*CZD[10,]^2-Vg_site)^2)/2)
Vg_crab <- (CZA[4,]*CZA[10,]^2+CZB[4,]*CZB[10,]^2+CZD[4,]*CZD[10,]^2)/3  # need to multiply by s^2 to get back to original scale
sd_Vgc <- sqrt(((CZA[4,]*CZA[10,]^2-Vg_crab)^2+(CZB[4,]*CZB[10,]^2-Vg_crab)^2+(CZD[4,]*CZD[10,]^2-Vg_crab)^2)/2)

params <- rbind(p, sd_p, Vg_site, sd_Vgs, Vg_crab, sd_Vgc) 

# get changes in phenotype in sd_g units
start <- means[1]
start_p <- start+p[phen]
change <- abs(means-start)/sqrt(as.numeric(Vg_crab[phen]))
change_p <- abs(means-as.numeric(start_p))/sqrt(as.numeric(Vg_crab[phen]))

print(data.frame(change,change_p))
