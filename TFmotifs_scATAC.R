##This script searches for your favorite TF motif in the "original" scATAC-seq differentially accessible peaks
##Established by Dana Bakalar
##Updated by Whitney Heavner 20240827

#Load the Genomic Ranges library
(library(GenomicRanges))

##Read in the rds file containing the motif positions. If you are using the .Rproj file system, this should be in your data directory.
df <- readRDS("data/Motif-Positions-In-Peaks.rds")

##Convert the rds file to a data frame
df2 <- as.data.frame(df)

##Read in the tables containing differentially accessible peaks in your cluster(s) of interest and change them to the Granges class
diff_acc <- as.data.frame(read.table("bed_files/Hybrid_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_range <- makeGRangesFromDataFrame(diff_acc)

diff_acc_act <- as.data.frame(read.table("bed_files/Activated_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_act_range <- makeGRangesFromDataFrame(diff_acc_act)

diff_acc_rest <- as.data.frame(read.table("bed_files/Resting_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_rest_range <- makeGRangesFromDataFrame(diff_acc_rest)

##Find ("subset") lines in the data frame containing your favorite TF motif
df_ebf <- df2[grep("Ebf",df2$group_name), ]
ebf_range <- makeGRangesFromDataFrame(df_ebf)

##Overlap the TF subset ranges with the ranges for differentially accessible peaks for your cluster(s) of interest
ebf_hybrid <- subsetByOverlaps(diff_acc_range, ebf_range)
ebf_act <- subsetByOverlaps(diff_acc_act_range, ebf_range)
ebf_rest <- subsetByOverlaps(diff_acc_rest_range, ebf_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
ebf_bed <- as.data.frame(ebf_hybrid)
ebf_bed <- ebf_bed[-c(5)]
write.table(ebf_bed, file="output/ebf_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

ebf_act_bed <- as.data.frame(ebf_act)
ebf_act_bed <- ebf_act_bed[-c(5)]
write.table(ebf_act_bed, file="output/ebf_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

ebf_rest_bed <- as.data.frame(ebf_rest)
ebf_rest_bed <- ebf_rest_bed[-c(5)]
write.table(ebf_rest_bed, file="output/ebf_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#Ebf containing peaks for hybrid cluster=479
#Ebf containing peaks for act cluster=1011
#Ebf containing peaks for rest cluster=311
#Total Ebf containing diff acc peaks=1801

#pval for hybrid cluster=phyper(total Ebf2 hybrid peaks, total Ebf2 peaks, total non-Ebf2 peaks, total hybrid peaks)
phyper(479, 1801, 24151-1801, 8492, lower.tail=FALSE, log.p=FALSE)
#1

#pval for activated cluster
phyper(1011, 1801, 24151-1801, 12146, lower.tail=FALSE, log.p=FALSE)
#1.070012e-07

#pval for resting cluster
phyper(311, 1801, 24151-1801, 3513, lower.tail=FALSE, log.p=FALSE)
#0.0003739926

##Foxa1 Example##

##Find ("subset") lines in the motif positions data frame containing the Foxa1 motif
df_fox <- df2[grep("Foxa1",df2$group_name), ]
foxa1_range <- makeGRangesFromDataFrame(df_fox)

##Overlap the TF subset with the differentially accessible peaks for the each cluster
foxa1_hybrid <- subsetByOverlaps(diff_acc_range, foxa1_range)
foxa1_act <- subsetByOverlaps(diff_acc_act_range, foxa1_range)
foxa1_rest <- subsetByOverlaps(diff_acc_rest_range, foxa1_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
foxa1_hyb_bed <- as.data.frame(foxa1_hybrid)
foxa1_hyb_bed <- foxa1_hyb_bed[-c(5)]
write.table(foxa1_hyb_bed, file="output/foxa1_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

foxa1_act_bed <- as.data.frame(foxa1_act)
foxa1_act_bed <- foxa1_act_bed[-c(5)]
write.table(foxa1_act_bed, file="output/foxa1_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

foxa1_rest_bed <- as.data.frame(foxa1_rest)
foxa1_rest_bed <- foxa1_rest_bed[-c(5)]
write.table(foxa1_rest_bed, file="output/foxa1_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

##This needs to be updated
#Foxa1 containing peaks for hybrid cluster=3061
#Foxa1 containing peaks for act cluster=969
#Foxa1 containing peaks for rest cluster=242
#Total Foxa1 containing diff acc peaks=4272

#pval for hybrid cluster=phyper(total Foxa1 hybrid peaks, total Foxa1 peaks, total non-Foxa1 peaks, total hybrid peaks)
phyper(3061, 4272, 24151-4272, 8492, lower.tail=FALSE, log.p=FALSE)
#0

#pval for activated cluster
phyper(969, 4272, 24151-4272, 12146, lower.tail=FALSE, log.p=FALSE)
#1

#pval for resting cluster
phyper(242, 4272, 24151-4272, 3513, lower.tail=FALSE, log.p=FALSE)
#1


##Ets Example##

##Find ("subset") lines in the motif positions data frame containing the motif
df_ets <- df2[grep("Ets",df2$group_name), ]
ets_range <- makeGRangesFromDataFrame(df_ets)

##Overlap the TF subset with the differentially accessible peaks for the each cluster
ets_hybrid <- subsetByOverlaps(diff_acc_range, ets_range)
ets_act <- subsetByOverlaps(diff_acc_act_range, ets_range)
ets_rest <- subsetByOverlaps(diff_acc_rest_range, ets_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
ets_hyb_bed <- as.data.frame(ets_hybrid)
ets_hyb_bed <- ets_hyb_bed[-c(5)]
write.table(ets_hyb_bed, file="output/ets_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

ets_act_bed <- as.data.frame(ets_act)
ets_act_bed <- ets_act_bed[-c(5)]
write.table(ets_act_bed, file="output/ets_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

ets_rest_bed <- as.data.frame(ets_rest)
ets_rest_bed <- ets_rest_bed[-c(5)]
write.table(ets_rest_bed, file="output/ets_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

##This needs to be updated
#Ets containing peaks for hybrid cluster=357
#Ets containing peaks for act cluster=582
#Ets containing peaks for rest cluster=134
#Total Ets containing diff acc peaks=1073

#pval for hybrid cluster=phyper(total ets hybrid peaks, total ets peaks, total non-ets peaks, total hybrid peaks)
phyper(357, 1073, 24151-1073, 8492, lower.tail=FALSE, log.p=FALSE)
#0.90

#pval for activated cluster
phyper(582, 1073, 24151-1073, 12146, lower.tail=FALSE, log.p=FALSE)
#0.0037

#pval for resting cluster
phyper(134, 1073, 24151-1073, 3513, lower.tail=FALSE, log.p=FALSE)
#0.97

##Klf5 Example##

##Find ("subset") lines in the motif positions data frame containing the motif
df_klf5 <- df2[grep("Klf5",df2$group_name), ]
klf5_range <- makeGRangesFromDataFrame(df_klf5)

##Overlap the TF subset with the differentially accessible peaks for the each cluster
klf5_hybrid <- subsetByOverlaps(diff_acc_range, klf5_range)
klf5_act <- subsetByOverlaps(diff_acc_act_range, klf5_range)
klf5_rest <- subsetByOverlaps(diff_acc_rest_range, klf5_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
klf5_hyb_bed <- as.data.frame(klf5_hybrid)
klf5_hyb_bed <- klf5_hyb_bed[-c(5)]
write.table(klf5_hyb_bed, file="output/klf5_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

klf5_act_bed <- as.data.frame(klf5_act)
klf5_act_bed <- klf5_act_bed[-c(5)]
write.table(klf5_act_bed, file="output/klf5_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

klf5_rest_bed <- as.data.frame(klf5_rest)
klf5_rest_bed <- klf5_rest_bed[-c(5)]
write.table(klf5_rest_bed, file="output/klf5_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#This needs to be updated
#klf5 containing peaks for hybrid cluster=1664
#klf5 containing peaks for act cluster=3114
#klf5 containing peaks for rest cluster=808
#Total klf5 containing diff acc peaks=5586

#pval for hybrid cluster=phyper(total klf5 hybrid peaks, total klf5 peaks, total non-klf5 peaks, total hybrid peaks)
phyper(1664, 5586, 24151-5586, 8492, lower.tail=FALSE, log.p=FALSE)
#1

#pval for activated cluster
phyper(3114, 5586, 24151-5586, 12146, lower.tail=FALSE, log.p=FALSE)
#5.6208e-21

#pval for resting cluster
phyper(808, 5586, 24151-5586, 3513, lower.tail=FALSE, log.p=FALSE)
#0.57

##Sox9 Example##

##Find ("subset") lines in the motif positions data frame containing the motif
df_sox9 <- df2[grep("Sox9",df2$group_name), ]
sox9_range <- makeGRangesFromDataFrame(df_sox9)

##Overlap the TF subset with the differentially accessible peaks for the each cluster
sox9_hybrid <- subsetByOverlaps(diff_acc_range, sox9_range)
sox9_act <- subsetByOverlaps(diff_acc_act_range, sox9_range)
sox9_rest <- subsetByOverlaps(diff_acc_rest_range, sox9_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
sox9_hyb_bed <- as.data.frame(sox9_hybrid)
sox9_hyb_bed <- sox9_hyb_bed[-c(5)]
write.table(sox9_hyb_bed, file="output/sox9_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

sox9_act_bed <- as.data.frame(sox9_act)
sox9_act_bed <- sox9_act_bed[-c(5)]
write.table(sox9_act_bed, file="output/sox9_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

sox9_rest_bed <- as.data.frame(sox9_rest)
sox9_rest_bed <- sox9_rest_bed[-c(5)]
write.table(sox9_rest_bed, file="output/sox9_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#sox9 containing peaks for hybrid cluster=519
#sox9 containing peaks for act cluster=629
#sox9 containing peaks for rest cluster=263
#Total sox9 containing diff acc peaks=1411

#pval for hybrid cluster=phyper(total sox9 hybrid peaks, total sox9 peaks, total non-sox9 peaks, total hybrid peaks)
phyper(519, 1411, 24151-1411, 8492, lower.tail=FALSE, log.p=FALSE)
#0.09

#pval for activated cluster
phyper(629, 1411, 24151-1411, 12146, lower.tail=FALSE, log.p=FALSE)
#0.99

#pval for resting cluster
phyper(263, 1411, 24151-1411, 3513, lower.tail=FALSE, log.p=FALSE)
#5.598962e-06


##Fos Example##
##Find ("subset") lines in the data frame containing your favorite TF motif
df_fos <- df2[grep("Fos",df2$group_name), ]
fos_range <- makeGRangesFromDataFrame(df_fos)

##Overlap the TF subset ranges with the ranges for differentially accessible peaks for your cluster(s) of interest
fos_hybrid <- subsetByOverlaps(diff_acc_range, fos_range)
fos_act <- subsetByOverlaps(diff_acc_act_range, fos_range)
fos_rest <- subsetByOverlaps(diff_acc_rest_range, fos_range)

##Write the new Granges list to a bed file; be sure to remove the 5th column first
fos_hyb_bed <- as.data.frame(fos_hybrid)
fos_hyb_bed <- fos_hyb_bed[-c(5)]
write.table(fos_hyb_bed, file="/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/fos_hybrid_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

fos_act_bed <- as.data.frame(fos_act)
fos_act_bed <- fos_act_bed[-c(5)]
write.table(fos_act_bed, file="/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/fos_act_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

fos_rest_bed <- as.data.frame(fos_rest)
fos_rest_bed <- fos_rest_bed[-c(5)]
write.table(fos_rest_bed, file="/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/fos_rest_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)

##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#Ebf containing peaks for hybrid cluster=479
#Ebf containing peaks for act cluster=1011
#Ebf containing peaks for rest cluster=311
#Total Ebf containing diff acc peaks=1801

#pval for hybrid cluster=phyper(total Ebf2 hybrid peaks, total Ebf2 peaks, total non-Ebf2 peaks, total hybrid peaks)
phyper(479, 1801, 24151-1801, 8492, lower.tail=FALSE, log.p=FALSE)
#1

#pval for activated cluster
phyper(1011, 1801, 24151-1801, 12146, lower.tail=FALSE, log.p=FALSE)
#1.070012e-07

#pval for resting cluster
phyper(311, 1801, 24151-1801, 3513, lower.tail=FALSE, log.p=FALSE)
#0.0003739926

#Overlap 2 Activated Peaks Bed Files (Fos and Ebf)##
#Find ("subset") lines in the data frame containing your favorite TF motif:
  #Fos
  df_fos <- df2[grep("Fos",df2$group_name), ]
  fos_range <- makeGRangesFromDataFrame(df_fos)
  
  #Ebf
  df_ebf <- df2[grep("Ebf",df2$group_name), ]
  ebf_range <- makeGRangesFromDataFrame(df_ebf)

#Overlap the TF subset ranges with the ranges for differentially accessible peaks for your cluster(s) of interest
  #Fos
  #fos_hybrid <- subsetByOverlaps(diff_acc_range, fos_range)
  fos_act <- subsetByOverlaps(diff_acc_act_range, fos_range)
  #fos_rest <- subsetByOverlaps(diff_acc_rest_range, fos_range)
  
  #Ebf
  #ebf_hybrid <- subsetByOverlaps(diff_acc_range, ebf_range)
  ebf_act <- subsetByOverlaps(diff_acc_act_range, ebf_range)
  #ebf_rest <- subsetByOverlaps(diff_acc_rest_range, ebf_range)
  
#alternatively, pull activated peaks data from previously generated bed files
  df_act_ebf <- as.data.frame(read.table("/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/ebf_act_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
  ebf_act <- makeGRangesFromDataFrame(df_act_ebf)
  
  df_act_fos <- as.data.frame(read.table("/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/fos_act_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
  ebf_act <- makeGRangesFromDataFrame(df_act_fos)
  
#Overlap the activated Fos peaks with activated Ebf peaks
fos_ebf_act <- subsetByOverlaps(fos_act, ebf_act)

#Write the new Granges list to a bed file; be sure to remove the 5th column first
fos_ebf_act_bed <- as.data.frame(fos_ebf_act)
fos_ebf_act_bed <- fos_ebf_act_bed[-c(5)]
write.table(fos_ebf_act_bed, file="/Volumes/MNSN/Ngai_Lab_Shared/10x/pipelines/TF_motif_finding/output/fos_ebf_activated_peaks.bed", quote=F, sep="\t", row.names=F, col.names=F)


##Calculate hypergeometic p-value for overlap for each cluster using lower tail=F
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#Ebf containing peaks for hybrid cluster=479
#Ebf containing peaks for act cluster=1011
#Ebf containing peaks for rest cluster=311
#Total Ebf containing diff acc peaks=1801

#pval for hybrid cluster=phyper(total Ebf2 hybrid peaks, total Ebf2 peaks, total non-Ebf2 peaks, total hybrid peaks)
phyper(479, 1801, 24151-1801, 8492, lower.tail=FALSE, log.p=FALSE)
#1

#pval for activated cluster
phyper(1011, 1801, 24151-1801, 12146, lower.tail=FALSE, log.p=FALSE)
#1.070012e-07

#pval for resting cluster
phyper(311, 1801, 24151-1801, 3513, lower.tail=FALSE, log.p=FALSE)
#0.0003739926

