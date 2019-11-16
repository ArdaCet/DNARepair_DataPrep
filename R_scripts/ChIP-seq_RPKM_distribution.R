library(tidyverse)
library(magrittr)
library(GGally)
library(ggplot2)
library(magrittr)
setwd('C:/Users/Arda/Documents/TEZ-Pairs-R')

chip_seqs <- read.table("99_new_ChIP-seqs_total.bed", header = FALSE, 
                        sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(chip_seqs)


### Selecting ChIP-seq data:
# "H3K79me2"	"H3K4me3"	"H3K4me3.1"	"H3K4me3.2"	"H3K27ac"	"H4K20me1"	"H3K27me3"	"H3K27me3.1"	"H3K36me3"	"H3K36me3.1"	"H3K4me2"	"H3K9ac"	"H3K9me3"	"H2AFZ"	"H3K4me1"

chip_seqs_no_control <- chip_seqs %>% filter(!V8 %in% c('Control', "rabbit-IgG-control", "mouse_IgG-control") )
head(chip_seqs_no_control)
unique(chip_seqs_no_control$V8)

chipseq_G1 = chip_seqs %>% filter(V8=="H3K79me2")
head(chipseq_G1)
chipseq_G2 = chip_seqs %>% filter(V8=="H3K4me3")

chipseq_G3 = chip_seqs %>% filter(V8=="H3K27ac") 

chipseq_G4 = chip_seqs %>% filter(V8=="H4K20me1") 

chipseq_G5 = chip_seqs %>% filter(V8=="H3K27me3") 

chipseq_G6 = chip_seqs %>% filter(V8=="H3K36me3") 

chipseq_G7 = chip_seqs %>% filter(V8=="H3K4me2") 

chipseq_G8 = chip_seqs %>% filter(V8=="H3K9ac") 

chipseq_G9 = chip_seqs %>% filter(V8=="H3K9me3") 

chipseq_G10 = chip_seqs %>% filter(V8=="H2AFZ") 

chipseq_G11 = chip_seqs %>% filter(V8=="H3K4me1")

# split G2 #

G2 = split(chipseq_G2, chipseq_G2$V9)
h3k4_1 <- G2[["SRR227441_SRR227442.fastq"]]
h3k4_2 <- G2[["SRR5338600_SRR5338596_SRR5338597_SRR5338598_SRR5338599.fastq"]] 
h3k4_3 <- G2[["SRR577378_SRR577379.fastq"]] 
head(h3k4_3)

# split G5 #
G5 = split(chipseq_G5, chipseq_G5$V9)
h3k27_1 <- G5[["SRR227473_SRR227472.fastq"]]
h3k27_2 <- G5[["SRR577392_SRR577393.fastq"]]

# split G6 #
G6 = split(chipseq_G6, chipseq_G6$V9)
h3k36_1 <- G6[["SRR227505_SRR227506.fastq"]]
h3k36_2 <- G6[["SRR577430_SRR577429.fastq"]]

rm(chipseq_G2, chipseq_G5, chipseq_G6)
rm(G2, G5, G6)
rm(chip_seqs)

###### TOTAL #########
c1 <- chipseq_G1 %>% mutate(H3K79me2 = V10) %>% select(H3K79me2)
c2_1 <- h3k4_1 %>% mutate(H3K4me3_1 = V10 ) %>% select(H3K4me3_1)
c2_2 <- h3k4_2 %>% mutate(H3K4me3_2 = V10 ) %>% select(H3K4me3_2)
c2_3 <- h3k4_3 %>% mutate(H3K4me3_3 = V10 ) %>% select(H3K4me3_3)
c3 <- chipseq_G3 %>% mutate( H3K27ac = V10 ) %>% select(H3K27ac)
c4 <- chipseq_G4 %>% mutate(H4K20me1 = V10 ) %>% select(H4K20me1)
c5_1 <- h3k27_1 %>% mutate(H3K27me3_1 = V10 ) %>% select(H3K27me3_1)
c5_2 <- h3k27_2 %>% mutate(H3K27me3_2 = V10 ) %>% select(H3K27me3_2)
c6_1 <- h3k36_1 %>% mutate(H3K36me3_1 = V10 ) %>% select(H3K36me3_1)
c6_2 <- h3k36_2 %>% mutate(H3K36me3_2 = V10 ) %>% select(H3K36me3_2)
c7 <- chipseq_G7 %>% mutate(H3K4me2 = V10 ) %>% select(H3K4me2)
c8 <- chipseq_G8 %>% mutate(H3K9ac = V10 ) %>% select(H3K9ac)
c9 <- chipseq_G9 %>% mutate(H3K9me3 = V10 ) %>% select(H3K9me3)
c10 <- chipseq_G10 %>% mutate(H2AFZ= V10) %>% select(H2AFZ)
c11 <- chipseq_G11 %>% mutate(H3K4me1 = V10) %>% select(H3K4me1)

rm(chipseq_G1,chipseq_G3, chipseq_G4, chipseq_G7,chipseq_G8, chipseq_G9,chipseq_G10,chipseq_G11,h3k4_1,h3k4_2,h3k4_3,h3k27_1,h3k27_2,h3k36_1,h3k36_2)

df_chip_seq = data.frame(c1, c2_1, c2_2, c2_3, c3, c4, c5_1, c5_2, c6_1, c6_2, c7, c8, c9, c10, c11)
head(df_chip_seq)

##################################  PLOTS ############################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15]) ### + 0.000001
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K79me2)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K79me2)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K79me2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_1)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K4me3_1)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_2)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K4me3_2)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_3)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K4me3_3)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_3_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27ac)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K27ac)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27ac_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H4K20me1)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H4K20me1)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H4K20me1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27me3_1)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K27me3_1)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27me3_2)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K27me3_2 )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K36me3_1 )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K36me3_1 )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K36me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K36me3_2 )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K36me3_2 )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K36me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me2 )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K4me2 )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K9ac )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K9ac )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K9ac_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K9me3 )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K9me3 )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K9me3_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H2AFZ )

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H2AFZ )) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H2AFZ_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15] )
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me1)

# Histogram Plots:
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Histogram")
### Specify jpeg height weight ###
w <- 595 
h <- 842
#Create Jpeg
p <- ggplot(df_chip_seq_log2, aes(x = H3K4me1)) + geom_histogram(bins=200) + scale_x_continuous(limits = c(0,30))
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)
rm(c1, c2_1, c2_2, c2_3, c3, c4, c5_1, c5_2, c6_1, c6_2, c7, c8, c9, c10, c11, df_chip_seq, chip_seqs_no_control)


#############################################################################################
############################### BOX PLOT  ###################################################
setwd("C:/Users/Arda/Documents/TEZ/TEZ-PlotResults/ChIP-seq_RPKM_distribution/log2_transformed/Box-plot")

# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K79me2)

p <- ggplot(df_chip_seq_log2, aes(y=H3K79me2, x="")) + stat_boxplot(aes("",H3K79me2) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K79me2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_1)

p <- ggplot(df_chip_seq_log2, aes(y=H3K4me3_1, x="")) + stat_boxplot(aes("",H3K4me3_1) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_2)

p <- ggplot(df_chip_seq_log2, aes(y=H3K4me3_2, x="")) + stat_boxplot(aes("",H3K4me3_2) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me3_3)

p <- ggplot(df_chip_seq_log2, aes(y=H3K4me3_3, x="")) + stat_boxplot(aes("",H3K4me3_3) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me3_3_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27ac)

p <- ggplot(df_chip_seq_log2, aes(y=H3K27ac, x="")) + stat_boxplot(aes("",H3K27ac) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27ac_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H4K20me1)

p <- ggplot(df_chip_seq_log2, aes(y=H4K20me1, x="")) + stat_boxplot(aes("",H4K20me1) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H4K20me1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27me3_1)

p <- ggplot(df_chip_seq_log2, aes(y=H3K27me3_1, x="")) + stat_boxplot(aes("",H3K27me3_1) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K27me3_2)

p <- ggplot(df_chip_seq_log2, aes(y=H3K27me3_2, x="")) + stat_boxplot(aes("",H3K27me3_2) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K27me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K36me3_1 )

p <- ggplot(df_chip_seq_log2, aes(y=H3K36me3_1, x="")) + stat_boxplot(aes("",H3K36me3_1) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K36me3_1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K36me3_2 )

p <- ggplot(df_chip_seq_log2, aes(y=H3K36me3_2, x="")) + stat_boxplot(aes("",H3K36me3_2) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K36me3_2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me2 )

p <- ggplot(df_chip_seq_log2, aes(y=H3K4me2, x="")) + stat_boxplot(aes("",H3K4me2) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me2_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K9ac )

p <- ggplot(df_chip_seq_log2, aes(y=H3K9ac, x="")) + stat_boxplot(aes("",H3K9ac) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K9ac_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K9me3 )

p <- ggplot(df_chip_seq_log2, aes(y=H3K9me3, x="")) + stat_boxplot(aes("",H3K9me3) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K9me3_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H2AFZ )

p <- ggplot(df_chip_seq_log2, aes(y=H2AFZ, x="")) + stat_boxplot(aes("",H2AFZ) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H2AFZ_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)

#############################################################################################
# log2 transformation:
df_chip_seq_log2 <- df_chip_seq

df_chip_seq_log2 <- log2(df_chip_seq_log2[1:15])
head(df_chip_seq_log2)
head(df_chip_seq_log2$H3K4me1)

p <- ggplot(df_chip_seq_log2, aes(y=H3K4me1, x="")) + stat_boxplot(aes("",H3K4me1) , geom='errorbar', linetype=1, width=0.05) + geom_boxplot()
p
#### Close jpeg ####
ggsave("log2transfromed_RPKMdistr_H3K4me1_noInfinity.jpg", p)
rm(p, df_chip_seq_log2)