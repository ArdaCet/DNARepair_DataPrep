library(tidyverse)
library(magrittr)
library(GGally)
library(ggplot2)

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

# Selecting columns for data fame:

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
                                                    H3K4me1_64PPB = as.numeric(V26)) %>% select(XR_HXA64B7_DMG, H3K79me2_64PPB, H3K4me3_64PPB, H3K4me3_1_64PPB, 
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
rm(RdD_1,RdD_2,RdD_3,RdD_4)                                                   

# Create Data-frames:

hmPPA = data.frame(RdD_1_filtered$XR_HXA64A1_DMG, RdD_1_filtered$H3K79me2_64PPA, RdD_1_filtered$H3K4me3_64PPA, RdD_1_filtered$H3K4me3_1_64PPA, RdD_1_filtered$H3K4me3_2_64PPA, RdD_1_filtered$H3K27ac_64PPA, RdD_1_filtered$H4K20me1_64PPA, RdD_1_filtered$H3K27me3_64PPA, RdD_1_filtered$H3K27me3_1_64PPA, RdD_1_filtered$H3K36me3_64PPA, RdD_1_filtered$H3K36me3_1_64PPA, RdD_1_filtered$H3K4me2_64PPA, RdD_1_filtered$H3K9ac_64PPA, RdD_1_filtered$H3K9me3_64PPA, RdD_1_filtered$H2AFZ_64PPA, RdD_1_filtered$H3K4me1_64PPA)
hmPPB = data.frame(RdD_2_filtered$XR_HXA64B7_DMG, RdD_2_filtered$H3K79me2_64PPB, RdD_2_filtered$H3K4me3_64PPB, RdD_2_filtered$H3K4me3_1_64PPB, RdD_2_filtered$H3K4me3_2_64PPB, RdD_2_filtered$H3K27ac_64PPB, RdD_2_filtered$H4K20me1_64PPB, RdD_2_filtered$H3K27me3_64PPB, RdD_2_filtered$H3K27me3_1_64PPB, RdD_2_filtered$H3K36me3_64PPB, RdD_2_filtered$H3K36me3_1_64PPB, RdD_2_filtered$H3K4me2_64PPB, RdD_2_filtered$H3K9ac_64PPB, RdD_2_filtered$H3K9me3_64PPB, RdD_2_filtered$H2AFZ_64PPB, RdD_2_filtered$H3K4me1_64PPB)
hmCPDA = data.frame(RdD_3_filtered$XR_HXACA4_DMG, RdD_3_filtered$H3K79me2_CPDA, RdD_3_filtered$H3K4me3_CPDA, RdD_3_filtered$H3K4me3_1_CPDA, RdD_3_filtered$H3K4me3_2_CPDA, RdD_3_filtered$H3K27ac_CPDA, RdD_3_filtered$H4K20me1_CPDA, RdD_3_filtered$H3K27me3_CPDA, RdD_3_filtered$H3K27me3_1_CPDA, RdD_3_filtered$H3K36me3_CPDA, RdD_3_filtered$H3K36me3_1_CPDA, RdD_3_filtered$H3K4me2_CPDA, RdD_3_filtered$H3K9ac_CPDA, RdD_3_filtered$H3K9me3_CPDA, RdD_3_filtered$H2AFZ_CPDA, RdD_3_filtered$H3K4me1_CPDA)
hmCPDB = data.frame(RdD_4_filtered$XR_HXACB10_DMG, RdD_4_filtered$H3K79me2_CPDB, RdD_4_filtered$H3K4me3_CPDB, RdD_4_filtered$H3K4me3_1_CPDB, RdD_4_filtered$H3K4me3_2_CPDB, RdD_4_filtered$H3K27ac_CPDB, RdD_4_filtered$H4K20me1_CPDB, RdD_4_filtered$H3K27me3_CPDB, RdD_4_filtered$H3K27me3_1_CPDB, RdD_4_filtered$H3K36me3_CPDB, RdD_4_filtered$H3K36me3_1_CPDB, RdD_4_filtered$H3K4me2_CPDB, RdD_4_filtered$H3K9ac_CPDB, RdD_4_filtered$H3K9me3_CPDB, RdD_4_filtered$H2AFZ_CPDB, RdD_4_filtered$H3K4me1_CPDB)

# Write down plots as PNG:
####################
setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K79me2_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K79me2_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K4me3_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K4me3_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K4me3_1_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K4me3_1_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K4me3_2_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K4me3_2_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K27ac_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K27ac_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H4K20me1_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H4K20me1_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K27me3_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K27me3_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K27me3_1_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K27me3_1_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K36me3_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K36me3_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K36me3_1_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K36me3_1_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K4me2_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K4me2_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K9ac_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K9ac_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K9me3_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K9me3_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H2AFZ_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H2AFZ_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64A1')
#Create Jpeg
png('H3K4me1_64PPA_wRegLine.png')
sp1 <- ggplot(hmPPA, aes(x=RdD_1_filtered$H3K4me1_64PPA, y=RdD_1_filtered$XR_HXA64A1_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)


############################################
setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')

#Create Jpeg
png('H3K79me2_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K79me2_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K4me3_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K4me3_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K4me3_1_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K4me3_1_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K4me3_2_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K4me3_2_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K27ac_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K27ac_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H4K20me1_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H4K20me1_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K27me3_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K27me3_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K27me3_1_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K27me3_1_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K36me3_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K36me3_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K36me3_1_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K36me3_1_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K4me2_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K4me2_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K9ac_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K9ac_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K9me3_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K9me3_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H2AFZ_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H2AFZ_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXA64B7')
#Create Jpeg
png('H3K4me1_64PPB_wRegLine.png')
sp1 <- ggplot(hmPPB, aes(x=RdD_2_filtered$H3K4me1_64PPB, y=RdD_2_filtered$XR_HXA64B7_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)


############################################
setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')

#Create Jpeg
png('H3K79me2_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K79me2_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K4me3_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K4me3_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K4me3_1_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K4me3_1_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K4me3_2_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K4me3_2_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K27ac_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K27ac_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H4K20me1_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H4K20me1_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K27me3_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K27me3_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K27me3_1_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K27me3_1_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K36me3_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K36me3_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K36me3_1_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K36me3_1_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K4me2_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K4me2_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K9ac_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K9ac_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K9me3_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K9me3_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H2AFZ_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H2AFZ_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACA4')
#Create Jpeg
png('H3K4me1_CPDA_wRegLine.png')
sp1 <- ggplot(hmCPDA, aes(x=RdD_3_filtered$H3K4me1_CPDA, y=RdD_3_filtered$XR_HXACA4_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)


############################################
setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')

#Create Jpeg
png('H3K79me2_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K79me2_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K4me3_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K4me3_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K4me3_1_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K4me3_1_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K4me3_2_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K4me3_2_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K27ac_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K27ac_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H4K20me1_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H4K20me1_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K27me3_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K27me3_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K27me3_1_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K27me3_1_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K36me3_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K36me3_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K36me3_1_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K36me3_1_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K4me2_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K4me2_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K9ac_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K9ac_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K9me3_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K9me3_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H2AFZ_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H2AFZ_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/Histone_x_Normalized-Repair-ScatterPlot/log-scaled-versions/HXACB10')
#Create Jpeg
png('H3K4me1_CPDB_wRegLine.png')
sp1 <- ggplot(hmCPDB, aes(x=RdD_4_filtered$H3K4me1_CPDB, y=RdD_4_filtered$XR_HXACB10_DMG)) + geom_point(color='blue', alpha = 0.01) + geom_smooth()
sp1 + scale_x_log10() + scale_y_log10()
#### Close jpeg ####
dev.off()
rm(sp1)
