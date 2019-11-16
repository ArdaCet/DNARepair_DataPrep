library(tidyverse)
library(magrittr)
library(GGally)

##################  xr_seq X Chip_seq ############################################

# import datasets as object type "data.frame"

xr_seqs <- read.table("C:/Users/Arda/Documents/TEZ-Pairs-R/total_xr-seqs.bed", header = FALSE, 
                   sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

#damage_seqs <- read.table("C:/Users/Arda/Documents/TEZ-Pairs-R/total_damage-seqs.bed", header = FALSE, 
#                      sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

chip_seqs <- read.table("C:/Users/Arda/Documents/TEZ-Pairs-R/99_new_ChIP-seqs_total.bed", header = FALSE, 
                        sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

# check if I got what I needed 

str(xr_seqs)
#str(damage_seqs)
str(chip_seqs)
unique(xr_seqs$V1) == unique(chip_seqs$V1)


# create the RPKM and Label columns for both datasets as featuring numeric and factor values respectively 

xr_seqs_filtered = xr_seqs %>% select(V1, V2, V3, V9, V10) %>% mutate(Start = V2, End = V3, 
                                                                      RPKM = as.numeric(V10), 
                                                                      Label = as.factor(V9), 
                                                                      CHR = as.factor(V1)) %>% select(CHR, Start, End, Label, RPKM)
#damage_seqs_filtered = damage_seqs %>% select(V1, V2, V3, V8, V9, V10, V11) %>% mutate(Start = V2, End = V3, 
#                                                                      RPKM = as.numeric(V11), 
#                                                                      Label = as.factor(V9),
#                                                                      Sign = as.factor(V10),
#                                                                      CHR = as.factor(V1)) %>% select(CHR, Start, End, Label, Sign, RPKM)
chip_seqs_filtered = chip_seqs %>% select(V1, V2, V3, V8, V9, V10) %>% mutate(Start = V2, End = V3, 
                                                                      RPKM = as.numeric(V10), 
                                                                      Group = as.factor(V8),
                                                                      Label = as.factor(V9),
                                                                      CHR = as.factor(V1)) %>% select(CHR, Start, End, Group, RPKM, Label)

# check again

str(xr_seqs_filtered)
#str(damage_seqs_filtered)
str(chip_seqs_filtered)
unique(chip_seqs_filtered$Group)
sum(chip_seqs_filtered$Start+5000 == chip_seqs_filtered$End)

# Make a column for the RPKM values of each group, make sure starts and ends match

chipseq_G1 = chip_seqs_filtered %>% filter(Group=="Control") %>% select(RPKM, Label) %>% mutate( Control = RPKM) %>% 
  select(Control, Label)

chipseq_G2 = chip_seqs_filtered %>% filter(Group=="H3K79me2") %>% select(RPKM) %>% mutate( H3K79me2 = RPKM) %>% 
  select(H3K79me2)

chipseq_G3 = chip_seqs_filtered %>% filter(Group=="H3K4me3") %>% select(RPKM, Label) %>% mutate( H3K4me3 = RPKM) %>% 
  select(H3K4me3, Label)

chipseq_G4 = chip_seqs_filtered %>% filter(Group=="H3K27ac") %>% select(RPKM) %>% mutate( H3K27ac = RPKM) %>% 
  select(H3K27ac)

chipseq_G5 = chip_seqs_filtered %>% filter(Group=="H4K20me1") %>% select(RPKM) %>% mutate( H4K20me1 = RPKM) %>% 
  select(H4K20me1)

chipseq_G6 = chip_seqs_filtered %>% filter(Group=="H3K27me3") %>% select(RPKM, Label) %>% mutate( H3K27me3 = RPKM) %>% 
  select(H3K27me3, Label)

chipseq_G7 = chip_seqs_filtered %>% filter(Group=="H3K36me3") %>% select(RPKM, Label) %>% mutate( H3K36me3 = RPKM) %>% 
  select(H3K36me3, Label)

chipseq_G8 = chip_seqs_filtered %>% filter(Group=="H3K4me2") %>% select(RPKM) %>% mutate( H3K4me2 = RPKM) %>% 
  select(H3K4me2)

chipseq_G9 = chip_seqs_filtered %>% filter(Group=="H3K9ac") %>% select(RPKM) %>% mutate( H3K9ac = RPKM) %>% 
  select(H3K9ac)

chipseq_G10 = chip_seqs_filtered %>% filter(Group=="rabbit-IgG-control") %>% select(RPKM) %>% 
  mutate(rabbit_IgG_control = RPKM) %>% select(rabbit_IgG_control)

chipseq_G11 = chip_seqs_filtered %>% filter(Group=="mouse_IgG-control") %>% select(RPKM) %>% 
  mutate( mouse_IgG_control = RPKM) %>% select(mouse_IgG_control)

chipseq_G12 = chip_seqs_filtered %>% filter(Group=="H3K9me3") %>% select(RPKM) %>% mutate( H3K9me3 = RPKM) %>% 
  select(H3K9me3)

chipseq_G13 = chip_seqs_filtered %>% filter(Group=="H2AFZ") %>% select(RPKM) %>% mutate( H2AFZ = RPKM) %>% 
  select(H2AFZ)

chipseq_G14 = chip_seqs_filtered %>% filter(Group=="H3K4me1") %>% select(RPKM) %>% mutate( H3K4me1 = RPKM) %>% 
  select(H3K4me1)


g1 = split(chipseq_G1, chipseq_G1$Label)
unique(chipseq_G1$Label)


### CONTROL GROUPS SPLIT ###
rm(chipseq_G1)
c1 <- g1[["SRR227392_SRR227391.fastq"]] %>% dplyr::select(Control)
c2 <- g1[["SRR351781.fastq"]] %>% dplyr::select(Control)
c3 <- g1[["SRR353505.fastq"]] %>% dplyr::select(Control)
c4 <- g1[["SRR353665.fastq"]] %>% dplyr::select(Control)
c5 <- g1[["SRR357521.fastq"]] %>% dplyr::select(Control)
c6 <- g1[["SRR502475_SRR502476.fastq"]] %>% dplyr::select(Control)
c7 <- g1[["SRR502477_SRR502478.fastq"]] %>% dplyr::select(Control)
c8 <- g1[["SRR5111806_SRR5111807.fastq"]] %>% dplyr::select(Control)
c9 <- g1[["SRR5112037_SRR5112036_SRR5112035.fastq"]] %>% dplyr::select(Control)
c10 <- g1[["SRR5331490_SRR5331491.fastq"]] %>% dplyr::select(Control)

class(g1[["SRR227392_SRR227391.fastq"]])

##### XR_SEQ SPLIT #######
head(xr_seqs_filtered)
unique(xr_seqs_filtered$Label)
xr_seqs_label1 = xr_seqs_filtered %>% filter(Label == "HXACA4_TGACCA") %>% select(RPKM) %>% mutate(HXACA4_TGACCA = RPKM) %>% select(HXACA4_TGACCA)
xr_seqs_label2 = xr_seqs_filtered %>% filter(Label == "HXACB10_TAGCTT") %>% select(RPKM) %>% mutate(HXACB10_TAGCTT = RPKM) %>% select(HXACB10_TAGCTT)
xr_seqs_label3 = xr_seqs_filtered %>% filter(Label == "HXA64B7_CAGATC") %>% select(RPKM) %>% mutate(HXA64B7_CAGATC = RPKM) %>% select(HXA64B7_CAGATC)
xr_seqs_label4 = xr_seqs_filtered %>% filter(Label == "HXA64A1_ATCACG") %>% select(RPKM) %>% mutate(HXA64A1_ATCACG = RPKM) %>% select(HXA64A1_ATCACG)

###### Damage_SEQ SPLIT ######
#head(damage_seqs_filtered)
#unique(damage_seqs_filtered$Label)
#damage_seqs_label1 = damage_seqs_filtered %>% filter(Label == "HDA64A1_ATCACG") %>% select(RPKM) %>% mutate(HDA64A1_ATCACG = RPKM) %>% select(HDA64A1_ATCACG)
#damage_seqs_label2 = damage_seqs_filtered %>% filter(Label == "HDA64B19_GTGAAA") %>% select(RPKM) %>% mutate(HDA64B19_GTGAAA = RPKM) %>% select(HDA64B19_GTGAAA)
#damage_seqs_label3 = damage_seqs_filtered %>% filter(Label == "HDACA6_GCCAAT") %>% select(RPKM) %>% mutate(HDACA6_GCCAAT = RPKM) %>% select(HDACA6_GCCAAT)
#damage_seqs_label4 = damage_seqs_filtered %>% filter(Label == "HDACB23_GAGTGG") %>% select(RPKM) %>% mutate(HDACB23_GAGTGG = RPKM) %>% select(HDACB23_GAGTGG)

###### G3 split ############
rm(g1)
g3 = split(chipseq_G3, chipseq_G3$Label)
str(chipseq_G3)

h3k4_1 <- g3[["SRR227441_SRR227442.fastq"]] %>% dplyr::select(H3K4me3)
h3k4_2 <- g3[["SRR5338600_SRR5338596_SRR5338597_SRR5338598_SRR5338599.fastq"]] %>% dplyr::select(H3K4me3)
h3k4_3 <- g3[["SRR577378_SRR577379.fastq"]] %>% dplyr::select(H3K4me3)

rm(g3)

##### G6 and G7 split ################
g6 = split(chipseq_G6, chipseq_G6$Label)
g7 = split(chipseq_G7, chipseq_G7$Label)
str(chipseq_G6)
str(chipseq_G7)
h3k27_1 <- g6[["SRR227473_SRR227472.fastq"]] %>% dplyr::select(H3K27me3)
h3k27_2 <- g6[["SRR577392_SRR577393.fastq"]] %>% dplyr::select(H3K27me3)

h3k36_1 <- g7[["SRR227505_SRR227506.fastq"]] %>% dplyr::select(H3K36me3)
h3k36_2 <- g7[["SRR577430_SRR577429.fastq"]] %>% dplyr::select(H3K36me3)

###### DELETE UNNECCESSARY STUFF #########
rm(g6)
rm(g7)
rm(chipseq_G6)
rm(chipseq_G7)
rm(chipseq_G3)
rm(xr_seqs)
rm(xr_seqs_filtered)
rm(chip_seqs)
rm(chip_seqs_filtered)
#rm(damage_seqs_filtered)
#rm(damage_seqs)

##### FINAL MERGE ###########
# c1, c2, c3, c4, c5, c6, c7, c8, c10 are removed from the data frame.
SIGMAR <- data.frame(chipseq_G2, chipseq_G4, chipseq_G5, chipseq_G8, chipseq_G9, 
           chipseq_G10, chipseq_G11, chipseq_G12, chipseq_G13, chipseq_G14, h3k27_1, h3k27_2, h3k36_1, h3k36_2,
           h3k4_1, h3k4_2, h3k4_3, xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4)

head(SIGMAR)
getwd()

round_1 = data.frame(xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4, chipseq_G2, chipseq_G4, chipseq_G5)
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

round_2 = data.frame(xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4, chipseq_G9, chipseq_G11, chipseq_G12)

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

round_3 = data.frame(xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4, chipseq_G13, chipseq_G14, h3k27_1, h3k27_2)

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

round_4 = data.frame(xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4, h3k4_1, h3k4_2, h3k4_3)

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

round_5 = data.frame(xr_seqs_label1, xr_seqs_label2, xr_seqs_label3, xr_seqs_label4, chipseq_G8 ,h3k36_1, h3k36_2)

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
jpeg('round5.jpg',width = w*10 , height = h*11)
pairs(round_5, log = "xy" , upper.panel=panel.cor, lower.panel = lower.panel)
#### Close jpeg ####
dev.off()

getwd()
