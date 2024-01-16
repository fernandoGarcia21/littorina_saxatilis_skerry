# Calculate allele frequency for each population & locus based on all data, as well as based on sub-
# sampling to the same sample size for all SNPs (needed for the demographic inference).
# Also write out arrangement frequencies for the inversions.
# The SNP input files contain the full raw data sets of control and spatial outlier SNPs used
# in this study. The columns give the individual genotypes, followed by
# cat = SNP ID for genotyping; contig = Contig; pos = Position in contig; LG = Linkage group;
# av = cM position in linkage group; cp = contig+position; inv_status = which inversion the SNP
# is located in (if any). The column IDs give the individual IDs; see below for how these
# translate to the different populations.
#
#author: "Anja Westram"
#
######################################################################################################


rm(list=ls())

# Set which loci to work with
type = "neu" # "neu" or "out" or "inv", for control SNPs / spatial outliers / inversions



# Get frequencies and frequencies after subsampling to same sample size for all loci
# for control and spatial SNPs 
if(type == "out" | type == "neu"){

  # Get genotype data
  if(type=="neu"){
    dat = read.table("../Data/SkerryExperiment_Neutral_NOLG12_Swedenfilter.txt", header=T, stringsAsFactors=F, sep=" ")}
  if(type=="out"){
    dat = read.table("../Data/SkerryExperiment_Outliers_NOLG12_Swedenfilter.txt", header=T, stringsAsFactors=F, sep=" ")}

  
  
  # Format genotype data frame
  rownames(dat) = dat$cp # Use locus ID as row name
  dat$cat = dat$contig = dat$pos = dat$LG = dat$av = dat$cp = dat$inv_status = NULL # Remove columns that are not needed
  dat[dat=="0/0"] = 0 # Reformat genotypes
  dat[dat=="0/1"] = 0.5
  dat[dat=="1/1"] = 1
  
  
  
  # Get population ID of each individual (C=Crab, W=Wave, S=Skerry; followed by year)
  pops = names(dat)
  pops[grep("DO.92", pops)] = "C1992"
  pops[grep("DonCrab", pops)] = "C2018"
  pops[grep("DO.21", pops)] = "C2021"
  pops[grep("RefW", pops)] = "W2018"
  pops[grep("W.REF.21", pops)] = "W2021"
  pops[grep("SKR.05", pops)] = "S2005"
  pops[grep("ExpSkerry", pops)] = "S2018"
  pops[grep("SRK.21", pops)] = "S2021"
  
  
  
  # Remove SNPs with too much missing data
  nas = apply(dat, 1, function(x) length(x[is.na(x)==T]))
  dat = dat[nas<50, ]
  
  
  
  # Make subsets for the different populations
  C1992 = dat[, pops=="C1992"]
  C0000 = dat[, pops %in% c("C1992", "C2018", "C2021")]
  W0000 = dat[, pops %in% c("W2018", "W2021")]
  S2005 = dat[, pops=="S2005"]
  S2018 = dat[, pops=="S2018"]
  S2021 = dat[, pops=="S2021"]
  
  
  
  # Get minimum number of genotypes per population (data will later be subsampled to this number)
  min(apply(C1992, 1, function(x) length(x[is.na(x)==F]))) # 18 19
  min(apply(C0000, 1, function(x) length(x[is.na(x)==F]))) #
  min(apply(W0000, 1, function(x) length(x[is.na(x)==F]))) # 83 77
  min(apply(S2005, 1, function(x) length(x[is.na(x)==F]))) # 31 32
  min(apply(S2018, 1, function(x) length(x[is.na(x)==F]))) # 45 48
  min(apply(S2021, 1, function(x) length(x[is.na(x)==F]))) # 36 34
  
  
  
  # Calculate allele frequencies
  pC1992 = apply(C1992, 1, function(x) mean(as.numeric(x), na.rm=T))
  pC0000 = apply(C0000, 1, function(x) mean(as.numeric(x), na.rm=T))
  pW0000 = apply(W0000, 1, function(x) mean(as.numeric(x), na.rm=T))
  pS2005 = apply(S2005, 1, function(x) mean(as.numeric(x), na.rm=T))
  pS2018 = apply(S2018, 1, function(x) mean(as.numeric(x), na.rm=T))
  pS2021 = apply(S2021, 1, function(x) mean(as.numeric(x), na.rm=T))
  
  
  
  # Make one data frame of all allele frequencies and remove fixed SNPs
  freqs = as.data.frame(cbind(dat$cp, pC1992, pC0000, pW0000, pS2005, pS2018, pS2021))
  keep = apply(freqs[, c("pC1992", "pC0000", "pW0000", "pS2005", "pS2018", "pS2021")],
               1, function(x) mean(as.numeric(x))!=0 & mean(as.numeric(x))!= 6)
  freqs = freqs[keep, ]
  
  
  
  # Write out
  freqs$cp = row.names(freqs)
  if(type=="neu"){
    write.table(freqs, "../Data/SEQSNPTM001_ALLELE_FREQS_NEU.txt", row.names=F, col.names=T, quote=F)}
  if(type=="out"){
    write.table(freqs, "../Data/SEQSNPTM001_ALLELE_FREQS_OUT.txt", row.names=F, col.names=T, quote=F)}
  
  
  
  # Now same steps again, but after downsampling to get equal sample size for each marker,
  # to allow for use with transition matrices (use minimum sample size across both spatial outlier
  # and control SNPs)
  pC1992 = apply(C1992, 1, function(x) mean(as.numeric(sample(x[is.na(x)==F], 18, replace=F))))
  pW0000 = apply(W0000, 1, function(x) mean(as.numeric(sample(x[is.na(x)==F], 77, replace=F))))
  pS2005 = apply(S2005, 1, function(x) mean(as.numeric(sample(x[is.na(x)==F], 31, replace=F))))
  pS2018 = apply(S2018, 1, function(x) mean(as.numeric(sample(x[is.na(x)==F], 45, replace=F))))
  pS2021 = apply(S2021, 1, function(x) mean(as.numeric(sample(x[is.na(x)==F], 34, replace=F))))
  
  freqs = as.data.frame(cbind(dat$cp, pC1992, pW0000, pS2005, pS2018, pS2021))
  
  keep = apply(freqs[, c("pC1992", "pW0000", "pS2005", "pS2018", "pS2021")],
               1, function(x) sum(as.numeric(x))!=0 & sum(as.numeric(x))!= 5)
  freqs = freqs[keep, ]
  
  freqs$cp = row.names(freqs)
  if(type=="neu"){
    write.table(freqs, "../Data/SEQSNPTM001_ALLELE_FREQS_NEU_subsampled.txt", row.names=F, col.names=T, quote=F)}
  if(type=="out"){
    write.table(freqs, "../Data/SEQSNPTM001_ALLELE_FREQS_OUT_subsampled.txt", row.names=F, col.names=T, quote=F)}
  
}



# Get frequencies for inversions and write out
# (This doesn't do anything except for changing one column name...)
if(type=="inv"){
  
  freqs2 = read.table("../Data/Inversions_Frequencies_for_Envelope_November2023.txt", header=T,
                      stringsAsFactors=F, sep="\t")
  names(freqs2)[1] = "cp"
  write.table(freqs2, "../Data/SEQSNPTM001_ALLELE_FREQS_INV.txt", row.names=F, col.names=T, quote=F)
  
}
