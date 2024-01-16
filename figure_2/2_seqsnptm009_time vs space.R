#
# author: "Anja Westram"
#
# Input file PAR001_freqs_ALL_13_all.txt is available on request
#
######################################################################################################

rm(list=ls())


# Function to calculate Fst
# Note: Use negative values where reference allele has lower frequency in Wave
# (This is done to be able to compare directionality of allele frequency change in space vs time)
fst = function(p1, p2){
  q1 = 1-p1
  q2 = 1-p2
  
  Ht = 2 * (p1+p2)/2 * (q1+q2)/2
  Hs = ((2*p1*q1) + (2*p2*q2)) / 2
  Fst = (Ht-Hs) / Ht
  Fst[is.na(Fst)] = 0
  Fst[p1>p2 & is.na(p1)==F & is.na(p2)==F] = -Fst[p1>p2 & is.na(p1)==F & is.na(p2)==F]
  
  return(Fst)
}


# Get skerry allele frequencies
obs = read.table("../Data/SEQSNPTM001_ALLELE_FREQS_OUT.txt",  header=T, stringsAsFactors=F)


# Get Fst over time on skerry
obs$fst_skerry = fst(obs$pC1992, obs$pS2021)


# Get outlier categories and add for each SNP
cat = read.table("../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt", header=T, stringsAsFactors=F)
cat$cat[grep("Hernan", cat$cat, value=F)] = "MoralesEtAl" # Outliers from Morales et al. 2019
cat$cat[grep("Anja", cat$cat, value=F)] = "WestramEtAl" # Outliers from Westram et al. 2018
obs = merge(obs, cat[, c("cp", "cat")], by="cp", all.x=T, all.y=F)


# Get pool-seq allele frequencies (Morales et al. 2019) and keep only nearby Swedish locations
# The file gives the allele frequencies in Crab (C) and Wave (W) populations from various locations;
# we keep only 4 nearby Swedish ones. Due to its large size, this file will be available only upon request.
pool = read.table("../Data/fst_previous_studies/PAR001_freqs_ALL_13_all.txt",
                 header=T, stringsAsFactors=F)
pool = pool[, c("cp","Sw_AC","Sw_AW","Sw_JC","Sw_JW","Sw_R3C","Sw_R3W","Sw_SC","Sw_SW")]


# Calculate Fst for each location and get average
pool = pool[is.na(pool$Sw_AC)==F, ]
pool$Sw_A = fst(pool$Sw_AC, pool$Sw_AW)
pool$Sw_J = fst(pool$Sw_JC, pool$Sw_JW)
pool$Sw_R3 = fst(pool$Sw_R3C, pool$Sw_R3W)
pool$Sw_S = fst(pool$Sw_SC, pool$Sw_SW)
pool$fst_area = apply(pool[,c("Sw_A","Sw_J","Sw_R3","Sw_S")], 1, function(x) mean(x, na.rm=T))


# Reverse all Fst values because the current study and the Morales et al. study used opposite reference alleles
pool$fst_area = -pool$fst_area


# Keep only outliers from Morales et al. 2019
MoralesEtAl = obs[obs$cat=="MoralesEtAl", ]


# Merge with poolseq data and orientate by spatial Fst (i.e. show allele with higher frequency in Wave)
MoralesEtAl = merge(MoralesEtAl, pool, all=F, by="cp")
rev = MoralesEtAl$fst_area<0
MoralesEtAl$fst_area[rev == T] = -MoralesEtAl$fst_area[rev == T] 
MoralesEtAl$fst_skerry[rev == T] = -MoralesEtAl$fst_skerry[rev == T] 


# Get data for 7 hybrid zones (Westram et al. 2021). In each hybrid zone file, each row is one SNP for
# which Fst between Crab and Wave was calculated. Only the column with the SNP ID ("cp") and the column
# with the Fst value ("FstRaw") are relevant
# Again, add sign to indicate directionality of differentiation
ang_right = read.table("../Data/fst_previous_studies/CZCLI006_ANG_rightFst.txt", header=T, stringsAsFactors=F)
ang_right$FstRaw[ang_right$p_crabRaw<ang_right$p_waveRaw] =
  - ang_right$FstRaw[ang_right$p_crabRaw<ang_right$p_waveRaw]
ang_right = ang_right[, c("cp","FstRaw")]

cza_left = read.table("../Data/fst_previous_studies/CZCLI006_CZA_leftFst.txt",
                       header=T, stringsAsFactors=F)
cza_left$FstRaw[cza_left$p_crabRaw<cza_left$p_waveRaw] =
  - cza_left$FstRaw[cza_left$p_crabRaw<cza_left$p_waveRaw]
cza_left = cza_left[, c("cp","FstRaw")]

cza_right = read.table("../Data/fst_previous_studies/CZCLI006_CZA_rightFst.txt",
                      header=T, stringsAsFactors=F)
cza_right$FstRaw[cza_right$p_crabRaw<cza_right$p_waveRaw] =
  - cza_right$FstRaw[cza_right$p_crabRaw<cza_right$p_waveRaw]
cza_right = cza_right[, c("cp","FstRaw")]

czb_left = read.table("../Data/fst_previous_studies/CZCLI006_CZB_leftFst.txt",
                      header=T, stringsAsFactors=F)
czb_left$FstRaw[czb_left$p_crabRaw<czb_left$p_waveRaw] =
  - czb_left$FstRaw[czb_left$p_crabRaw<czb_left$p_waveRaw]
czb_left = czb_left[, c("cp","FstRaw")]

czb_right = read.table("../Data/fst_previous_studies/CZCLI006_CZB_rightFst.txt",
                       header=T, stringsAsFactors=F)
czb_right$FstRaw[czb_right$p_crabRaw<czb_right$p_waveRaw] =
  - czb_right$FstRaw[czb_right$p_crabRaw<czb_right$p_waveRaw]
czb_right = czb_right[, c("cp","FstRaw")]

czd_left = read.table("../Data/fst_previous_studies/CZCLI006_CZD_leftFst.txt",
                      header=T, stringsAsFactors=F)
czd_left$FstRaw[czd_left$p_crabRaw<czd_left$p_waveRaw] =
  - czd_left$FstRaw[czd_left$p_crabRaw<czd_left$p_waveRaw]
czd_left = czd_left[, c("cp","FstRaw")]

czd_right = read.table("../Data/fst_previous_studies/CZCLI006_CZD_rightFst.txt",
                       header=T, stringsAsFactors=F)
czd_right$FstRaw[czd_right$p_crabRaw<czd_right$p_waveRaw] =
  - czd_right$FstRaw[czd_right$p_crabRaw<czd_right$p_waveRaw]
czd_right = czd_right[, c("cp","FstRaw")]


# Get average Fst in space
cz = as.data.frame(cbind(ang_right, cza_left$FstRaw, cza_right$FstRaw, czb_left$FstRaw,
           czb_right$FstRaw, czd_left$FstRaw, czd_right$FstRaw))
cz$fst_area = apply(cz[, 2:8], 1, function(x) mean(x, na.rm=T))


# Keep only outliers from Westram et al. 2018
WestramEtAl = obs[obs$cat=="WestramEtAl", ]


# Merge with cz data and orientate by spatial Fst (i.e. show allele with higher frequency in Wave)
WestramEtAl = merge(WestramEtAl, cz, all=F, by="cp")
rev = WestramEtAl$fst_area<0
WestramEtAl$fst_area[rev == T] = -WestramEtAl$fst_area[rev == T] 
WestramEtAl$fst_skerry[rev == T] = -WestramEtAl$fst_skerry[rev == T] 


# Plot space-time correlation
MoralesEtAl$pch = 16
plot(MoralesEtAl$fst_area, MoralesEtAl$fst_skerry, col=rgb(0.2,0.6,0.6,0.4), pch=MoralesEtAl$pch,
     xlab="Crab-Wave Fst in space (nearby locations)", ylab="Fst in time (Skerry)")

WestramEtAl$pch = 16
points(WestramEtAl$fst_area, WestramEtAl$fst_skerry, col=rgb(0,1,1,0.4), pch=WestramEtAl$pch)


# Get model and add line to plot
mod = lm(c(WestramEtAl$fst_skerry, MoralesEtAl$fst_skerry) ~ c(WestramEtAl$fst_area, MoralesEtAl$fst_area))
abline(mod, lwd=2)
summary(mod)
cor.test(c(WestramEtAl$fst_skerry, MoralesEtAl$fst_skerry), c(WestramEtAl$fst_area, MoralesEtAl$fst_area),
         method = "spearman")


# Write out
res = rbind(MoralesEtAl[, c("cp", "cat", "fst_area","fst_skerry")],
            WestramEtAl[, c("cp", "cat", "fst_area","fst_skerry")])
write.table(res, "../Data/SEQSNPTM009_TIME VS SPACE_OUT.txt", col.names=T, row.names=F, quote=F)




# Use similar approach for inversions

# Get skerry dataset
invs = read.table(paste("../Data/SEQSNPTM001_ALLELE_FREQS_", "INV", ".txt", sep=""),
                 header=T, stringsAsFactors=F)


# Keep only arrangement most common in Wave for the complex inversions
invs = invs[duplicated(invs$cp)==F, ]


# Get temporal Fst
invs$fst_skerry = fst(invs$pC1992, invs$pS2021)
     

# Get contact zone dataset, which shows arrangement frequencies in the 7 locations from Westram et al. 2021
czinvs = read.table("../Data/Table_inversion_frequencies_Diego_final_simplified.csv",
                    header=T, sep=",", stringsAsFactors=F)


# The data are orientated to show the allele most common in the Wave reference population, so change to
# the other allele where needed (i.e. where the skerry data shows the other allele)
invs$cp[invs$pW0000<0.5]
czinvs[czinvs$Inversions=="LGC9.1", 2:15] = 1-czinvs[czinvs$Inversions=="LGC9.1", 2:15]
czinvs[czinvs$Inversions=="LGC17.1", 2:15] = 1-czinvs[czinvs$Inversions=="LGC17.1", 2:15]


# Get average Fst across locations
for(zone in c("ANG", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")){
  czinvs[, zone] = fst(czinvs[, paste(zone, "CRAB", sep="_")], czinvs[, paste(zone, "Wave", sep="_")])
}
czinvs$fst_area = apply(czinvs[, c("ANG", "CZA_left", "CZA_right", "CZB_left", "CZB_right", "CZD_left", "CZD_right")],
                        1, function(x) mean(x, na.rm=T))


# Combine spatial and temporal Fsts
invs = merge(invs, czinvs[, c("Inversions", "fst_area")], by.x="cp", by.y="Inversions")


# Orientate by spatial Fst
rev = invs$fst_area<0
invs$fst_area[rev == T] = -invs$fst_area[rev == T] 
invs$fst_skerry[rev == T] = -invs$fst_skerry[rev == T] 


# Plot
plot(invs$fst_area, invs$fst_skerry, col="red", pch=16)
text(invs$fst_area, invs$fst_skerry, labels=invs$cp, pos=1, cex=0.5) 


# Get model and add line to plot
mod = lm(invs$fst_skerry ~ invs$fst_area)
abline(mod, lwd=2)
summary(mod)
cor.test(invs$fst_skerry, invs$fst_area,
         method = "spearman")


# Write out
write.table(invs[, c("cp", "fst_skerry", "fst_area")], "../Data/SEQSNPTM009_TIME VS SPACE_INV.txt",
            col.names=T, row.names=F, quote=F)