library(tidyverse)
library(magrittr)
library(GGally)

# Read Files:
setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/XR-seqs_ChIP-seqs-divided-noZero-Data')

all_seq1 <- read.table("HXA64A1_ATCACG_ChiP-seqs_tabdelimited_noZeroRows_divided_DF1.txt",skip = 1, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

all_seq2 <- read.table("HXA64B7_CAGATC_ChiP-seqs_tabdelimited_noZeroRows_divided_DF2.txt",skip = 1, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

all_seq3 <- read.table("HXACA4_TGACCA_ChiP-seqs_tabdelimited_noZeroRows_divided_DF3.txt",skip = 1, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

all_seq4 <- read.table("HXACB10_TAGCTT_ChiP-seqs_tabdelimited_noZeroRows_divided_DF4.txt",skip = 1, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(all_seq1)
head(all_seq2)
head(all_seq3)
head(all_seq4)

#Filter Files:
# DF1
xr_seqs_label1 = all_seq1 %>% select(V1) %>% mutate(HXA64A1_ATCACG = V1) %>% select(HXA64A1_ATCACG)
c1 <- all_seq1 %>% mutate(H3K79me2 = V11) %>% select(H3K79me2)
c2_1 <- all_seq1 %>% mutate(H3K4me3 = V21 ) %>% select(H3K4me3)
c2_2 <- all_seq1 %>% mutate(H3K4me3 = V22 ) %>% select(H3K4me3)
c2_3 <- all_seq1 %>% mutate(H3K4me3 = V23 ) %>% select(H3K4me3)
c3 <- all_seq1 %>% mutate( H3K27ac = V12 ) %>% select(H3K27ac)
c4 <- all_seq1 %>% mutate(H4K20me1 = V13 ) %>% select(H4K20me1)
c5_1 <- all_seq1 %>% mutate(H3K27me3 = V24 ) %>% select(H3K27me3)
c5_2 <- all_seq1 %>% mutate(H3K27me3 = V25 ) %>% select(H3K27me3)
c6_1 <- all_seq1 %>% mutate(H3K36me3 = V26 ) %>% select(H3K36me3)
c6_2 <- all_seq1 %>% mutate(H3K36me3 = V27 ) %>% select(H3K36me3)
c7 <- all_seq1 %>% mutate(H3K4me2 = V14 ) %>% select(H3K4me2)
c8 <- all_seq1 %>% mutate(H3K9ac = V15 ) %>% select(H3K9ac)
c9 <- all_seq1 %>% mutate(H3K9me3 = V18 ) %>% select(H3K9me3)
c10 <- all_seq1 %>% mutate(H2AFZ= V19) %>% select(H2AFZ)
c11 <- all_seq1 %>% mutate(H3K4me1 = V20) %>% select(H3K4me1)
# DF2
xr_seqs_label2 = all_seq2 %>% select(V1) %>% mutate(HXA64B7_CAGATC = V1) %>% select(HXA64B7_CAGATC)
d1 <- all_seq2 %>% mutate(H3K79me2 = V11) %>% select(H3K79me2)
d2_1 <- all_seq2 %>% mutate(H3K4me3 = V21 ) %>% select(H3K4me3)
d2_2 <- all_seq2 %>% mutate(H3K4me3 = V22 ) %>% select(H3K4me3)
d2_3 <- all_seq2 %>% mutate(H3K4me3 = V23 ) %>% select(H3K4me3)
d3 <- all_seq2 %>% mutate( H3K27ac = V12 ) %>% select(H3K27ac)
d4 <- all_seq2 %>% mutate(H4K20me1 = V13 ) %>% select(H4K20me1)
d5_1 <- all_seq2 %>% mutate(H3K27me3 = V24 ) %>% select(H3K27me3)
d5_2 <- all_seq2 %>% mutate(H3K27me3 = V25 ) %>% select(H3K27me3)
d6_1 <- all_seq2 %>% mutate(H3K36me3 = V26 ) %>% select(H3K36me3)
d6_2 <- all_seq2 %>% mutate(H3K36me3 = V27 ) %>% select(H3K36me3)
d7 <- all_seq2 %>% mutate(H3K4me2 = V14 ) %>% select(H3K4me2)
d8 <- all_seq2 %>% mutate(H3K9ac = V15 ) %>% select(H3K9ac)
d9 <- all_seq2 %>% mutate(H3K9me3 = V18 ) %>% select(H3K9me3)
d10 <- all_seq2 %>% mutate(H2AFZ= V19) %>% select(H2AFZ)
d11 <- all_seq2 %>% mutate(H3K4me1 = V20) %>% select(H3K4me1)
# DF3
xr_seqs_label3 = all_seq3 %>% select(V1) %>% mutate(HXACA4_TGACCA = V1) %>% select(HXACA4_TGACCA)
e1 <- all_seq3 %>% mutate(H3K79me2 = V11) %>% select(H3K79me2)
e2_1 <- all_seq3 %>% mutate(H3K4me3 = V21 ) %>% select(H3K4me3)
e2_2 <- all_seq3 %>% mutate(H3K4me3 = V22 ) %>% select(H3K4me3)
e2_3 <- all_seq3 %>% mutate(H3K4me3 = V23 ) %>% select(H3K4me3)
e3 <- all_seq3 %>% mutate( H3K27ac = V12 ) %>% select(H3K27ac)
e4 <- all_seq3 %>% mutate(H4K20me1 = V13 ) %>% select(H4K20me1)
e5_1 <- all_seq3 %>% mutate(H3K27me3 = V24 ) %>% select(H3K27me3)
e5_2 <- all_seq3 %>% mutate(H3K27me3 = V25 ) %>% select(H3K27me3)
e6_1 <- all_seq3 %>% mutate(H3K36me3 = V26 ) %>% select(H3K36me3)
e6_2 <- all_seq3 %>% mutate(H3K36me3 = V27 ) %>% select(H3K36me3)
e7 <- all_seq3 %>% mutate(H3K4me2 = V14 ) %>% select(H3K4me2)
e8 <- all_seq3 %>% mutate(H3K9ac = V15 ) %>% select(H3K9ac)
e9 <- all_seq3 %>% mutate(H3K9me3 = V18 ) %>% select(H3K9me3)
e10 <- all_seq3 %>% mutate(H2AFZ= V19) %>% select(H2AFZ)
e11 <- all_seq3 %>% mutate(H3K4me1 = V20) %>% select(H3K4me1)
# DF4
xr_seqs_label4 = all_seq4 %>% select(V1) %>% mutate(HXACB10_TAGCTT = V1) %>% select(HXACB10_TAGCTT)
f1 <- all_seq4 %>% mutate(H3K79me2 = V11) %>% select(H3K79me2)
f2_1 <- all_seq4 %>% mutate(H3K4me3 = V21 ) %>% select(H3K4me3)
f2_2 <- all_seq4 %>% mutate(H3K4me3 = V22 ) %>% select(H3K4me3)
f2_3 <- all_seq4 %>% mutate(H3K4me3 = V23 ) %>% select(H3K4me3)
f3 <- all_seq4 %>% mutate( H3K27ac = V12 ) %>% select(H3K27ac)
f4 <- all_seq4 %>% mutate(H4K20me1 = V13 ) %>% select(H4K20me1)
f5_1 <- all_seq4 %>% mutate(H3K27me3 = V24 ) %>% select(H3K27me3)
f5_2 <- all_seq4 %>% mutate(H3K27me3 = V25 ) %>% select(H3K27me3)
f6_1 <- all_seq4 %>% mutate(H3K36me3 = V26 ) %>% select(H3K36me3)
f6_2 <- all_seq4 %>% mutate(H3K36me3 = V27 ) %>% select(H3K36me3)
f7 <- all_seq4 %>% mutate(H3K4me2 = V14 ) %>% select(H3K4me2)
f8 <- all_seq4 %>% mutate(H3K9ac = V15 ) %>% select(H3K9ac)
f9 <- all_seq4 %>% mutate(H3K9me3 = V18 ) %>% select(H3K9me3)
f10 <- all_seq4 %>% mutate(H2AFZ= V19) %>% select(H2AFZ)
f11 <- all_seq4 %>% mutate(H3K4me1 = V20) %>% select(H3K4me1)
#############################################################################
rm(all_seq1,all_seq2,all_seq3,all_seq4)
# Start Pairs Ploting:

# For DF1:
round_1 = data.frame(xr_seqs_label1,c1,c2_1,c2_2,c2_3)
head(round_1)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))  # Pearson_default
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round1.jpg',width = w*10 , height = h*11)
pairs(round_1, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_2 = data.frame(xr_seqs_label1,c3,c4,c5_1,c5_2)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round2.jpg',width = w*10 , height = h*11)
pairs(round_2, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_3 = data.frame(xr_seqs_label1,c6_1,c6_2,c7,c8)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round3.jpg',width = w*10 , height = h*11)
pairs(round_3, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_4 = data.frame(xr_seqs_label1,c9,c10,c11)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round4.jpg',width = w*10 , height = h*11)
pairs(round_4, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

#################################################################################
# For DF2:

round_6 = data.frame(xr_seqs_label2,d1,d2_1,d2_2,d2_3)
head(round_6)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))  # Pearson_default
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round6.jpg',width = w*10 , height = h*11)
pairs(round_6, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_7 = data.frame(xr_seqs_label2,d3,d4,d5_1,d5_2)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round7.jpg',width = w*10 , height = h*11)
pairs(round_7, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_8 = data.frame(xr_seqs_label2,d6_1,d6_2,d7,d8)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round8.jpg',width = w*10 , height = h*11)
pairs(round_8, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_9 = data.frame(xr_seqs_label2,d9,d10,d11)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round9.jpg',width = w*10 , height = h*11)
pairs(round_9, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

############################################################################
# For DF3:

round_11 = data.frame(xr_seqs_label3,e1,e2_1,e2_2,e2_3)
head(round_11)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))  # Pearson_default
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round11.jpg',width = w*10 , height = h*11)
pairs(round_11, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_12 = data.frame(xr_seqs_label3,e3,e4,e5_1,e5_2)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round12.jpg',width = w*10 , height = h*11)
pairs(round_12, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_13 = data.frame(xr_seqs_label3,e6_1,e6_2,e7,e8)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round13.jpg',width = w*10 , height = h*11)
pairs(round_13, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_14 = data.frame(xr_seqs_label3,e9,e10,e11)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round14.jpg',width = w*10 , height = h*11)
pairs(round_14, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

#########################################################################
# For DF4:

round_16 = data.frame(xr_seqs_label4,f1,f2_1,f2_2,f2_3)
head(round_16)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))  # Pearson_default
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round16.jpg',width = w*10 , height = h*11)
pairs(round_16, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_17 = data.frame(xr_seqs_label4,f3,f4,f5_1,f5_2)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round17.jpg',width = w*10 , height = h*11)
pairs(round_17, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_18 = data.frame(xr_seqs_label4,f6_1,f6_2,f7,f8)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round18.jpg',width = w*10 , height = h*11)
pairs(round_18, log = "xy" ,upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

round_19 = data.frame(xr_seqs_label4,f9,f10,f11)

color=rgb(0,0.2,0.8,alpha=0.3) 

# Correlation panel
panel.cor <- panel.cor <- function(x, y, ...)
{
  par(usr = c(0, 1, 0, 1), xlog = FALSE, ylog = FALSE)
  txt <- as.character(format(cor(x, y), digits=2))
  text(0.5, 0.5, txt, cex = 6* abs(cor(x, y)))
}
# Customize lower panel
lower.panel<-function(x, y){
  points(x,y, pch = 16, col = color)
}
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
jpeg('round19.jpg',width = w*10 , height = h*11)
pairs(round_19, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

##########################################################################