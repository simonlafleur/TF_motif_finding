##This script searches for your favorite TF motif in the "original" scATAC-seq differentially accessible peaks
##Established by Dana Bakalar
##Updated by Whitney Heavner 20240827

#Load the Genomic Ranges library
(library(GenomicRanges))

##Read in the rds file containing the motif positions. If you are using the .Rproj file system, this should be in your data directory.
df <- readRDS("data/Motif-Positions-In-Peaks.rds")

##Convert the rds file to a data frame
df2 <- as.data.frame(df)

##Find ("subset") lines in the data frame containing your favorite TF motif
df_ebf <- df2[grep("Ebf",df2$group_name), ]
ebf_range <- makeGRangesFromDataFrame(df_ebf)

##Read in the tables containing differentially accessible peaks in your cluster(s) of interest and change them to the Granges class
diff_acc <- as.data.frame(read.table("bed_files/Hybrid_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_range <- makeGRangesFromDataFrame(diff_acc)

diff_acc_act <- as.data.frame(read.table("bed_files/Activated_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_act_range <- makeGRangesFromDataFrame(diff_acc_act)

diff_acc_rest <- as.data.frame(read.table("bed_files/Resting_diff_acc_peaks.bed", header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
diff_acc_rest_range <- makeGRangesFromDataFrame(diff_acc_rest)

##Overlap the TF subset with the differentially accessible peaks for your cluster(s) of interest
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

##Calculate hypergeometic p-value for overlap for each cluster
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#Ebf containing peaks for hybrid cluster=543
#Ebf containing peaks for act cluster=1029
#Ebf containing peaks for rest cluster=317
#Total Ebf containing diff acc peaks=1889

#pval for hybrid cluster=phyper(total Ebf2 hybrid peaks, total Ebf2 peaks, total non-Ebf2 peaks, total hybrid peaks)
phyper(543, 1889, 24151-1889, 8492, lower.tail=FALSE, log.p=FALSE)
#1

#pval for activated cluster
phyper(1029, 1889, 24151-1889, 12146, lower.tail=FALSE, log.p=FALSE)
#6.886279e-05

#pval for resting cluster
phyper(317, 1889, 24151-1889, 3513, lower.tail=FALSE, log.p=FALSE)
#0.002136152

##Foxa1 Example##

##Find ("subset") lines in the motif positions data frame containing the Uncx motif
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

##Calculate hypergeometic p-value for overlap for each cluster
#Diff acc peaks total for hybrid cluster=8492
#Diff acc peaks total for act cluster=12146
#Diff acc peaks total for rest cluster=3513
#Total diff acc peaks=24151

#Foxa1 containing peaks for hybrid cluster=3265
#Foxa1 containing peaks for act cluster=1020
#Foxa1 containing peaks for rest cluster=251
#Total Foxa1 containing diff acc peaks=4536

#pval for hybrid cluster=phyper(total Foxa1 hybrid peaks, total Foxa1 peaks, total non-Foxa1 peaks, total hybrid peaks)
phyper(3265, 4536, 24151-4536, 8492, lower.tail=FALSE, log.p=FALSE)
#0

#pval for activated cluster
phyper(1020, 4536, 24151-4536, 12146, lower.tail=FALSE, log.p=FALSE)
#1

#pval for resting cluster
phyper(251, 4536, 24151-4536, 3513, lower.tail=FALSE, log.p=FALSE)
#1

