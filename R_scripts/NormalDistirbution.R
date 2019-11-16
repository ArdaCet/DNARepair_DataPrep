library(ggplot2)
library(tidyverse)
library(magrittr)
library(GGally)

                                                   #  For stock DMG files:

setwd("C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/Damage-seq_Data")

DF1 <- read.table('HDA64A1_ATCACG_DF1_RPKMsumation-added.csv', skip = 2,header = FALSE, sep="\t",
                     stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(DF1)

DF2 <- read.table('HDA64B19_GTGAAA_DF2_RPKMsumation-added.csv', skip = 2, header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(DF2)

DF3 <- read.table('HDACA6_GCCAAT_DF3_RPKMsumation-added.csv', skip = 2, header = FALSE, sep="\t",
                 stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(DF3)

DF4 = read.table('HDACB23_GAGTGG_DF4_RPKMsumation-added.csv', skip = 2, header = FALSE, sep="\t",
                 stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/NormalDistributions')

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF1:
p <- ggplot(DF1, aes(x=V1)) + labs(y="Number of RPKM values", x="6-4PP replicate A") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDA64A1_ATCACG_histogram_DF1.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF2:
p <- ggplot(DF2, aes(x=V1)) + labs(y="Number of RPKM values", x="6-4PP replicate B") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDA64B19_GTGAAA_histogram_DF2.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF3:
p <- ggplot(DF3, aes(x=V1)) + labs(y="Number of RPKM values", x="CPD replicate A") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDACA6_GCCAAT_histogram_DF3.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF4:
p <- ggplot(DF4, aes(x=V1)) + labs(y="Number of RPKM values", x="CPD replicate B") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDACB23_GAGTGG_histogram_DF4.jpg', p)
rm(p)

## BACKUP:
#############################################
qplot(x = HDACA6_GCCAAT_M, data= y3_m, 
      bins = 600,
      geom="histogram", 
      main = "Histogram for Damage-seq", 
      xlab = "RPKM", 
      fill=I("blue"))

qplot(x = HDACA6_GCCAAT_P, data= y3_p, 
      bins = 600,
      geom="histogram", 
      main = "Histogram for Damage-seq", 
      xlab = "RPKM", 
      fill=I("red"))

###---------------------------------------------------------------------------------------------------------------------------
rm(DF1, DF2, DF3, DF4)
                                                #  For Zero removed DMG files:


setwd("C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/Damage-seq_Data/PlusandMinus_added_DMG_WithoutZeros")

DF1_n0 <- read.table('HDA64A1_ATCACG_noZeroRows_new_DF1.csv', header = FALSE, sep="\t",
           stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(DF1_n0)
DF2_n0 = read.table('HDA64B19_GTGAAA_noZeroRows_new_DF2.csv', header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

DF3_n0 = read.table('HDACA6_GCCAAT_noZeroRows_new_DF3.csv', header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
   
DF4_n0 = read.table('HDACB23_GAGTGG_noZeroRows_new_DF4.csv', header = FALSE, sep="\t",
                    stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)


setwd('C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/NormalDistributions')

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF1:
p <- ggplot(DF1_n0, aes(x=V1)) + labs(y="Number of RPKM values", x="6-4PP replicate A") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDA64A1_ATCACG_histogram_noZero_DF1.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF2:
p <- ggplot(DF2_n0, aes(x=V1)) + labs(y="Number of RPKM values", x="6-4PP replicate B") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDA64B19_GTGAAA_histogram_noZero_DF2.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF3:
p <- ggplot(DF3_n0, aes(x=V1)) + labs(y="Number of RPKM values", x="CPD replicate A") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDACA6_GCCAAT_histogram_noZero_DF3.jpg', p)
rm(p)

### Specify jpeg height weight ###
w <- 595 
h <- 842
# Plot DF4:
p <- ggplot(DF4_n0, aes(x=V1)) + labs(y="Number of RPKM values", x="CPD replicate B") + geom_histogram(binwidth = 10) + xlim(-5,7000)
p
#### Close jpeg ####
ggsave('HDACB23_GAGTGG_histogram_noZero_DF4.jpg', p)
rm(p)
