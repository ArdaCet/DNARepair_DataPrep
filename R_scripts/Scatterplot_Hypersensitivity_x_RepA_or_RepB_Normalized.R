library(tidyverse)
library(magrittr)
library(GGally)
library(ggplot2)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/Peak-Repair-DMG-ChIP-seq_wo_Cntrl-DATA/Only_RPKM_Data/STOCK-Data_wo_Cntrl/WithCHRStartEnd')

# ChIP-seq File without Controls:

ChIP_fl <- read.table("XR-DMG-ChIP_wo_Cntrl-Peak_Data_w_ChrStartEnd.csv", header = FALSE, 
                   sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(ChIP_fl)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/Peak-Repair-DMG-ChIP-seq_wo_Cntrl-DATA/Only_RPKM_Data/NoZeroRow-divided/All_seq_Peak_divided_zeroRemoved-DATA/wChr-Start-End')

# 6-4PP Damage's XR-seq, Damage-seq and ChIP-seq Data:

RdD_1 <- read.table("ZeroRmd_XR_over_DMG_all-seq-including_Peak_HXA64A1_w_ChrStartEnd.csv", header = FALSE, 
                    sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_1)
RdD_2 <- read.table("ZeroRmd_XR_over_DMG_all-seq-including_Peak_HXA64B19_w_ChrStartEnd.csv", header = FALSE, 
                    sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_2)

# CPD Damage's XR-seq, Damage-seq and ChIP-seq Data:

RdD_3 <- read.table("ZeroRmd_XR_over_DMG_all-seq-including_Peak_HXACA6_w_ChrStartEnd.csv", header = FALSE, 
                    sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_3)
RdD_4 <- read.table("ZeroRmd_XR_over_DMG_all-seq-including_Peak_HXACB23_w_ChrStartEnd.csv", header = FALSE, 
                    sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_4)


setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/ReplicateA-Normalized-woZero-Divided-allseq-Data')

# 6-4PP Damage's REP_A Normalized XR-seq, Damage-seq and ChIP-seq Data:

RdD_5 <- read.csv("Zero_removed_divided_HDA64A1_ATCACG_normalized_allseqs_w_ChrStartEnd.csv", header = FALSE, 
                  sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_5)

# CPD Damage's REP_A Normalized  XR-seq, Damage-seq and ChIP-seq Data:

RdD_6 <- read.csv("Zero_removed_divided_HDACA6_GCCAAT_normalized_allseqs_w_ChrStartEnd.csv", header = FALSE, 
                  sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_6)


setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/ReplicateB-Normalized-woZero-Divided-allseq-Data')

# 6-4PP Damage's REP_B Normalized XR-seq, Damage-seq and ChIP-seq Data:

RdD_7 <- read.csv("Zero_removed_divided_HDA64B19_GTGAAA_normalized_allseqs_w_ChrStartEnd.csv", header = FALSE, 
                  sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_7)

# CPD Damage's REP_B Normalized XR-seq, Damage-seq and ChIP-seq Data:

RdD_8 <- read.csv("Zero_removed_divided_HDACB23_GAGTGG_normalized_allseqs_w_ChrStartEnd.csv", header = FALSE, 
                  sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(RdD_8)


setwd("D:/SRR/DNase-I-Hyersensitivity-regions/bedfile")

# DNase I Hyperaensitivity Sites, DNase-seq Datum:

DIhS <- read.table("99_SRR352412_SRR352413_SRR352414.fastq_RPKM_added.bed", header = FALSE, 
                   sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(DIhS)


setwd("C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/DNaseI_Hypersensitivity_Sites/Zero_rows/Normalized_Zero-row")

# DNase I Hypersensitivity Sites filtered, DNase-seq Datum:

DIhs_PPD_RepA <- read.table("Zero_removed_HDA64A1_ATCACG_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv", header = FALSE, 
                            sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(DIhs_PPD_RepA)

DIhs_PPD_RepB <- read.table("Zero_removed_HDA64B19_GTGAAA_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv", header = FALSE, 
                            sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(DIhs_PPD_RepB)

DIhs_CPD_RepA <- read.table("Zero_removed_HDACA6_GCCAAT_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv", header = FALSE, 
                            sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(DIhs_CPD_RepA)

DIhs_CPD_RepB <- read.table("Zero_removed_HDACB23_GAGTGG_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv", header = FALSE, 
                            sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(DIhs_CPD_RepB)
# Selecting columns for data fame:
ChIP_fl_filtered = ChIP_fl %>% select(V12:V26) %>% mutate(H3K79me2 = as.numeric(V12), H3K4me3 = as.numeric(V13), 
                                                          H3K4me3_1 = as.numeric(V14), H3K4me3_2 = as.numeric(V15), H3K27ac = as.numeric(V16), 
                                                          H4K20me1 = as.numeric(V17), H3K27me3 = as.numeric(V18), H3K27me3_1 = as.numeric(V19), 
                                                          H3K36me3 = as.numeric(V20), H3K36me3_1 = as.numeric(V21), H3K4me2 = as.numeric(V22), 
                                                          H3K9ac = as.numeric(V23), H3K9me3 = as.numeric(V24), H2AFZ = as.numeric(V25), 
                                                          H3K4me1 = as.numeric(V26)) %>% select(H3K79me2, H3K4me3, H3K4me3_1, 
                                                                                                      H3K4me3_2, H3K27ac, H4K20me1, 
                                                                                                      H3K27me3, H3K27me3_1, H3K36me3, 
                                                                                                      H3K36me3_1, H3K4me2, H3K9ac, H3K9me3, 
                                                                                                      H2AFZ, H3K4me1)


RdD_1_filtered = RdD_1 %>% select(V1:V26) %>% mutate(XR_HXA64A1_DMG = as.numeric(V4), H3K79me2_64PPA = as.numeric(V12), H3K4me3_64PPA = as.numeric(V13), 
                                                     H3K4me3_1_64PPA = as.numeric(V14), H3K4me3_2_64PPA = as.numeric(V15), H3K27ac_64PPA = as.numeric(V16), 
                                                     H4K20me1_64PPA = as.numeric(V17), H3K27me3_64PPA = as.numeric(V18), H3K27me3_1_64PPA = as.numeric(V19), 
                                                     H3K36me3_64PPA = as.numeric(V20), H3K36me3_1_64PPA = as.numeric(V21), H3K4me2_64PPA = as.numeric(V22), 
                                                     H3K9ac_64PPA = as.numeric(V23), H3K9me3_64PPA = as.numeric(V24), H2AFZ_64PPA = as.numeric(V25), 
                                                     H3K4me1_64PPA = as.numeric(V26)) %>% select(XR_HXA64A1_DMG, H3K79me2_64PPA, H3K4me3_64PPA, H3K4me3_1_64PPA, 
                                                                                                 H3K4me3_2_64PPA, H3K27ac_64PPA, H4K20me1_64PPA, 
                                                                                                 H3K27me3_64PPA, H3K27me3_1_64PPA, H3K36me3_64PPA, 
                                                                                                 H3K36me3_1_64PPA, H3K4me2_64PPA, H3K9ac_64PPA, H3K9me3_64PPA, 
                                                                                                 H2AFZ_64PPA, H3K4me1_64PPA)

RdD_2_filtered = RdD_2 %>% select(V1:V26) %>% mutate(XR_HXA64B7_DMG = as.numeric(V5), H3K79me2_64PPB = as.numeric(V12), H3K4me3_64PPB = as.numeric(V13), 
                                                     H3K4me3_1_64PPB = as.numeric(V14), H3K4me3_2_64PPB = as.numeric(V15), H3K27ac_64PPB = as.numeric(V16), 
                                                     H4K20me1_64PPB = as.numeric(V17), H3K27me3_64PPB = as.numeric(V18), H3K27me3_1_64PPB = as.numeric(V19), 
                                                     H3K36me3_64PPB = as.numeric(V20), H3K36me3_1_64PPB = as.numeric(V21), H3K4me2_64PPB = as.numeric(V22), 
                                                     H3K9ac_64PPB = as.numeric(V23), H3K9me3_64PPB = as.numeric(V24), H2AFZ_64PPB = as.numeric(V25), 
                                                     H3K4me1_64PPB = as.numeric(V26)) %>% select(XR_HXA64B7_DMG,H3K79me2_64PPB, H3K4me3_64PPB, H3K4me3_1_64PPB, 
                                                                                                 H3K4me3_2_64PPB, H3K27ac_64PPB, H4K20me1_64PPB, 
                                                                                                 H3K27me3_64PPB, H3K27me3_1_64PPB, H3K36me3_64PPB, 
                                                                                                 H3K36me3_1_64PPB, H3K4me2_64PPB, H3K9ac_64PPB, H3K9me3_64PPB,
                                                                                                 H2AFZ_64PPB, H3K4me1_64PPB)

RdD_3_filtered = RdD_3 %>% select(V1:V26) %>% mutate(XR_HXACA4_DMG = as.numeric(V6), H3K79me2_CPDA = as.numeric(V12), H3K4me3_CPDA = as.numeric(V13), 
                                                     H3K4me3_1_CPDA = as.numeric(V14), H3K4me3_2_CPDA = as.numeric(V15), H3K27ac_CPDA = as.numeric(V16), 
                                                     H4K20me1_CPDA = as.numeric(V17), H3K27me3_CPDA = as.numeric(V18), H3K27me3_1_CPDA = as.numeric(V19), 
                                                     H3K36me3_CPDA = as.numeric(V20), H3K36me3_1_CPDA = as.numeric(V21), H3K4me2_CPDA = as.numeric(V22), 
                                                     H3K9ac_CPDA = as.numeric(V23), H3K9me3_CPDA = as.numeric(V24), H2AFZ_CPDA = as.numeric(V25), 
                                                     H3K4me1_CPDA = as.numeric(V26)) %>% select(XR_HXACA4_DMG, H3K79me2_CPDA, H3K4me3_CPDA, H3K4me3_1_CPDA, 
                                                                                                H3K4me3_2_CPDA, H3K27ac_CPDA, H4K20me1_CPDA, 
                                                                                                H3K27me3_CPDA, H3K27me3_1_CPDA, H3K36me3_CPDA, 
                                                                                                H3K36me3_1_CPDA, H3K4me2_CPDA, H3K9ac_CPDA, H3K9me3_CPDA, 
                                                                                                H2AFZ_CPDA, H3K4me1_CPDA)

RdD_4_filtered = RdD_4 %>% select(V1:V26) %>% mutate(XR_HXACB10_DMG = as.numeric(V7), H3K79me2_CPDB = as.numeric(V12), H3K4me3_CPDB = as.numeric(V13), 
                                                     H3K4me3_1_CPDB = as.numeric(V14), H3K4me3_2_CPDB = as.numeric(V15), H3K27ac_CPDB = as.numeric(V16), 
                                                     H4K20me1_CPDB = as.numeric(V17), H3K27me3_CPDB = as.numeric(V18), H3K27me3_1_CPDB = as.numeric(V19), 
                                                     H3K36me3_CPDB = as.numeric(V20), H3K36me3_1_CPDB = as.numeric(V21), H3K4me2_CPDB = as.numeric(V22), 
                                                     H3K9ac_CPDB = as.numeric(V23), H3K9me3_CPDB = as.numeric(V24), H2AFZ_CPDB = as.numeric(V25), 
                                                     H3K4me1_CPDB = as.numeric(V26)) %>% select(XR_HXACB10_DMG, H3K79me2_CPDB, H3K4me3_CPDB, H3K4me3_1_CPDB, 
                                                                                                H3K4me3_2_CPDB, H3K27ac_CPDB, H4K20me1_CPDB, 
                                                                                                H3K27me3_CPDB, H3K27me3_1_CPDB, H3K36me3_CPDB, 
                                                                                                H3K36me3_1_CPDB, H3K4me2_CPDB, H3K9ac_CPDB, H3K9me3_CPDB, 
                                                                                                H2AFZ_CPDB, H3K4me1_CPDB)

RdD_5_filtered = RdD_5 %>% select(V4) %>% mutate(XR_HXA64A1_RepA_DMG = as.numeric(V4)) %>% select(XR_HXA64A1_RepA_DMG)
RdD_6_filtered = RdD_5 %>% select(V5) %>% mutate(XR_HXA64B7_RepA_DMG = as.numeric(V5)) %>% select(XR_HXA64B7_RepA_DMG)
RdD_7_filtered = RdD_6 %>% select(V4) %>% mutate(XR_HXACA4_RepA_DMG = as.numeric(V4)) %>% select(XR_HXACA4_RepA_DMG)
RdD_8_filtered = RdD_6 %>% select(V5) %>% mutate(XR_HXACB10_RepA_DMG = as.numeric(V5)) %>% select(XR_HXACB10_RepA_DMG)

RdD_9_filtered = RdD_7 %>% select(V4) %>% mutate(XR_HXA64A1_RepB_DMG = as.numeric(V4)) %>% select(XR_HXA64A1_RepB_DMG) 
RdD_10_filtered = RdD_7 %>% select(V5) %>% mutate(XR_HXA64B7_RepB_DMG = as.numeric(V5)) %>% select(XR_HXA64B7_RepB_DMG)
RdD_11_filtered = RdD_8 %>% select(V4) %>% mutate(XR_HXACA4_RepB_DMG = as.numeric(V4)) %>% select(XR_HXACA4_RepB_DMG)
RdD_12_filtered = RdD_8 %>% select(V5) %>% mutate(XR_HXACB10_RepB_DMG = as.numeric(V5)) %>% select(XR_HXACB10_RepB_DMG)

DIhS_filtered = DIhS %>% select(V10) %>% mutate(DNase_I_Hypersensitivity = as.numeric(V10)) %>% select(DNase_I_Hypersensitivity)

DIhs_PPD_RepA_filtered = DIhs_PPD_RepA %>% select(V10) %>% mutate(DNase_I_Hypersensitivity_PPD_RepA = as.numeric(V10)) %>% select(DNase_I_Hypersensitivity_PPD_RepA)
DIhs_PPD_RepB_filtered = DIhs_PPD_RepB %>% select(V10) %>% mutate(DNase_I_Hypersensitivity_PPD_RepB = as.numeric(V10)) %>% select(DNase_I_Hypersensitivity_PPD_RepB)
DIhs_CPD_RepA_filtered = DIhs_CPD_RepA %>% select(V10) %>% mutate(DNase_I_Hypersensitivity_CPD_RepA = as.numeric(V10)) %>% select(DNase_I_Hypersensitivity_CPD_RepA)
DIhs_CPD_RepB_filtered = DIhs_CPD_RepB %>% select(V10) %>% mutate(DNase_I_Hypersensitivity_CPD_RepB = as.numeric(V10)) %>% select(DNase_I_Hypersensitivity_CPD_RepB)

rm(RdD_1,RdD_2,RdD_3,RdD_4,RdD_5,RdD_6,RdD_7,RdD_8,DIhS,DIhs_PPD_RepA,DIhs_PPD_RepB,DIhs_CPD_RepA,DIhs_CPD_RepB,ChIP_fl)                                                   

# Create Data-frames:
# Normalized Filtered Repair and DNase I hypersensitivity regions
hmPPA = data.frame(RdD_1_filtered$XR_HXA64A1_DMG, DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)
hmPPB = data.frame(RdD_2_filtered$XR_HXA64B7_DMG, DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)
hmCPDA = data.frame(RdD_3_filtered$XR_HXACA4_DMG, DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)
hmCPDB = data.frame(RdD_4_filtered$XR_HXACB10_DMG, DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)

# Replicate A Normalized Repair and DNase I hypersensitivity regions
hmPPA_repA = data.frame(RdD_5_filtered$XR_HXA64A1_RepA_DMG, DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)
hmPPB_repA = data.frame(RdD_6_filtered$XR_HXA64B7_RepA_DMG, DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)
hmCPDA_repA = data.frame(RdD_7_filtered$XR_HXACA4_RepA_DMG, DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)
hmCPDB_repA = data.frame(RdD_8_filtered$XR_HXACB10_RepA_DMG, DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)

# Replicate B Normalized Repair and DNase I hypersensitivity regions
hmPPA_repB = data.frame(RdD_9_filtered$XR_HXA64A1_RepB_DMG, DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)
hmPPB_repB = data.frame(RdD_10_filtered$XR_HXA64B7_RepB_DMG,DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)
hmCPDA_repB = data.frame(RdD_11_filtered$XR_HXACA4_RepB_DMG,DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)
hmCPDB_repB = data.frame(RdD_12_filtered$XR_HXACB10_RepB_DMG,DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)

# Histone marker and DNase I hypersensitivity regions
hmDIh = data.frame(DIhS_filtered$DNase_I_Hypersensitivity,ChIP_fl_filtered$H3K79me2, ChIP_fl_filtered$H3K4me3, ChIP_fl_filtered$H3K4me3_1, ChIP_fl_filtered$H3K4me3_2, ChIP_fl_filtered$H3K27ac, ChIP_fl_filtered$H4K20me1, ChIP_fl_filtered$H3K27me3, ChIP_fl_filtered$H3K27me3_1, ChIP_fl_filtered$H3K36me3, ChIP_fl_filtered$H3K36me3_1, ChIP_fl_filtered$H3K4me2, ChIP_fl_filtered$H3K9ac, ChIP_fl_filtered$H3K9me3, ChIP_fl_filtered$H2AFZ, ChIP_fl_filtered$H3K4me1)

########################################################################################################################




# Write down plots as PNG:

##  FOR HISTONE MARKERS and DNase I HYPERSENSITIVITY SCATTERPLOT
####################
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K79me2_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K79me2, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K4me3_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K4me3, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K4me3_1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K4me3_1, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K4me3_2_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K4me3_2, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K27ac_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K27ac, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H4K20me1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H4K20me1, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K27me3_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K27me3, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K27me3_1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K27me3_1, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K36me3_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K36me3, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K36me3_1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K36me3_1, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K4me2_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K4me2, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K9ac_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K9ac, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K9me3_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K9me3, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H2AFZ_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H2AFZ, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Histone_x_DNaseI")
#Create Jpeg
png('H3K4me1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmDIh, aes(x=ChIP_fl_filtered$H3K4me1, y=DIhS_filtered$DNase_I_Hypersensitivity)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)


##  FOR FILTERED and NORMALIZED XR-SEQ SCATTERPLOT
####################
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Norm_Repair_x_DNaseI")
#Create Jpeg
png('HXA64A1_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$XR_HXA64A1_DMG, y=DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Norm_Repair_x_DNaseI")

#Create Jpeg
png('HXA64B7_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$XR_HXA64B7_DMG, y=DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Norm_Repair_x_DNaseI")

#Create Jpeg
png('HXACA4_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$XR_HXACA4_DMG, y=DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Norm_Repair_x_DNaseI")

#Create Jpeg
png('HXACB10_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$XR_HXACB10_DMG, y=DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

##  FOR FILTERED and REPLICATE A NORMALIZED XR-SEQ SCATTERPLOT
####################
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_A_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXA64A1_RepA_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPA_repA, aes(x=RdD_5_filtered$XR_HXA64A1_RepA_DMG, y=DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_A_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXA64B7_RepA_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPB_repA, aes(x=RdD_6_filtered$XR_HXA64B7_RepA_DMG, y=DIhs_PPD_RepA_filtered$DNase_I_Hypersensitivity_PPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_A_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXACA4_RepA_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDA_repA, aes(x=RdD_7_filtered$XR_HXACA4_RepA_DMG, y=DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_A_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXACB10_RepA_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDB_repA, aes(x=RdD_8_filtered$XR_HXACB10_RepA_DMG, y=DIhs_CPD_RepA_filtered$DNase_I_Hypersensitivity_CPD_RepA)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

##  FOR FILTERED and REPLICATE B NORMALIZED XR-SEQ SCATTERPLOT
####################
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_B_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXA64A1_RepB_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPA_repB, aes(x=RdD_9_filtered$XR_HXA64A1_RepB_DMG, y=DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_B_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXA64B7_RepB_DNaseI_wRegLine.png')
sp1 <- ggplot(hmPPB_repB, aes(x=RdD_10_filtered$XR_HXA64B7_RepB_DMG, y=DIhs_PPD_RepB_filtered$DNase_I_Hypersensitivity_PPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_B_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXACA4_RepB_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDA_repB, aes(x=RdD_11_filtered$XR_HXACA4_RepB_DMG, y=DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_DNaseI-Hypersensitivity-ScatterPlot/Rep_B_Norm_Repair_x_DNaseI")

#Create Jpeg
png('XR_HXACB10_RepB_DNaseI_wRegLine.png')
sp1 <- ggplot(hmCPDB_repB, aes(x=RdD_12_filtered$XR_HXACB10_RepB_DMG, y=DIhs_CPD_RepB_filtered$DNase_I_Hypersensitivity_CPD_RepB)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)
