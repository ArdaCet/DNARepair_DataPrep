library(tidyverse)
library(magrittr)
library(GGally)

# XR-seq and ChIP-seq Data
setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Stock_Seq_Files/Subsampled')

xr_seq <- read.table("subsampled_GM_xrseq.bed", header = FALSE, 
                    sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(xr_seq)

chip_seq <- read.table("subsampled_GM_chipseq.bed", header = FALSE, 
                       sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(chip_seq)

#DNaseI <- read.table("DNaseI-seq_total.bed", header = FALSE, 
#                     sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
#head(DNaseI)

# Damage-seq Data wo Sign diversity
setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Damage-seq_Data/No_ControlandSignDiversity_Data/Subsampled')

dmg_seq <- read.table("dmg_seqs_no_signs_subsampled_woHeader.csv", header = FALSE, 
                      sep="\t",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(dmg_seq)

##############################################################################################################################
# Filter ChIP-seq Data HeLa cells:

unique(chip_seq$V8)

h4k20 = chip_seq %>% filter(V8=="H4K20me1") %>% select(V1,V2,V3,V10) %>% mutate(H4K20me1_name = V1, H4K20me1_start = V2, H4K20me1_end = V3, H4K20me1_RPKM = V10) %>% select(H4K20me1_RPKM)

h3k9me3 = chip_seq %>% filter(V8=="H3K9me3") %>% select(V1,V2,V3,V10) %>% mutate(H3K9me3_name = V1, H3K9me3_start = V2, H3K9me3_end = V3, H3K9me3_RPKM = V10) %>% select(H3K9me3_RPKM)

h3k9ac  = chip_seq %>% filter(V8=="H3K9ac") %>% select(V1,V2,V3,V10) %>% mutate(H3K9ac_name = V1, H3K9ac_start = V2, H3K9ac_end = V3, H3K9ac_RPKM = V10) %>% select(H3K9ac_RPKM)

h3k79me2  = chip_seq %>% filter(V8=="H3K79me2") %>% select(V1,V2,V3,V10) %>% mutate(H3K79me2_name = V1, H3K79me2_start = V2, H3K79me2_end = V3, H3K79me2_RPKM = V10) %>% select(H3K79me2_RPKM)

h2afz  = chip_seq %>% filter(V8=="H2AFZ") %>% select(V1,V2,V3,V10) %>% mutate(H2AFZ_name = V1, H2AFZ_start = V2, H2AFZ_end = V3, H2AFZ_RPKM = V10) %>% select(H2AFZ_RPKM)

h3k4me2  = chip_seq %>% filter(V8=="H3K4me2") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me2_name = V1, H3K4me2_start = V2, H3K4me2_end = V3, H3K4me2_RPKM = V10) %>% select(H3K4me2_RPKM)

h3k27ac  = chip_seq %>% filter(V8=="H3K27ac") %>% select(V1,V2,V3,V10) %>% mutate(H3K27ac_name = V1, H3K27ac_start = V2, H3K27ac_end = V3, H3K27ac_RPKM = V10) %>% select(H3K27ac_RPKM)

h3k4me1  = chip_seq %>% filter(V8=="H3K4me1") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me1_name = V1, H3K4me1_start = V2, H3K4me1_end = V3, H3K4me1_RPKM = V10) %>% select(H3K4me1_RPKM)

h3k36me3  = chip_seq %>% filter(V9=="SRR577430_SRR577429.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K36me3_name = V1, H3K36me3_start = V2, H3K36me3_end = V3, H3K36me3_RPKM = V10) %>% select(H3K36me3_RPKM)

h3k36me3_2  = chip_seq %>% filter(V9=="SRR227505_SRR227506.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K36me3_2_name = V1, H3K36me3_2_start = V2, H3K36me3_2_end = V3, H3K36me3_2_RPKM = V10) %>% select(H3K36me3_2_RPKM)

h3k27me3  = chip_seq %>% filter(V9=="SRR577392_SRR577393.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K27me3_name = V1, H3K27me3_start = V2, H3K27me3_end = V3, H3K27me3_RPKM = V10) %>% select(H3K27me3_RPKM)

h3k27me3_2  = chip_seq %>% filter(V9=="SRR227473_SRR227472.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K27me3_2_name = V1, H3K27me3_2_start = V2, H3K27me3_2_end = V3, H3K27me3_2_RPKM = V10) %>% select(H3K27me3_2_RPKM)

h3k4me3  = chip_seq %>% filter(V9=="SRR227441_SRR227442.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me3_name = V1, H3K4me3_start = V2, H3K4me3_end = V3, H3K4me3_RPKM = V10) %>% select(H3K4me3_RPKM)

h3k4me3_2  = chip_seq %>% filter(V9=="SRR5338600_SRR5338596_SRR5338597_SRR5338598_SRR5338599.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me3_2_name = V1, H3K4me3_2_start = V2, H3K4me3_2_end = V3, H3K4me3_2_RPKM = V10) %>% select(H3K4me3_2_RPKM)

h3k4me3_3 = chip_seq %>% filter(V9 =="SRR577378_SRR577379.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me3_3_name = V1, H3K4me3_3_start = V2, H3K4me3_3_end = V3, H3K4me3_3_RPKM = V10) %>% select(H3K4me3_3_RPKM)

############################################################################################################################
# Filter ChIP-seq Data GM12878 cells:

unique(chip_seq$V8)

h3k4me3  = chip_seq %>% filter(V9=="SRR5331171_SRR5331172.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me3_name = V1, H3K4me3_start = V2, H3K4me3_end = V3, H3K4me3_RPKM = V10) %>% select(H3K4me3_RPKM)

h3k4me3_2  = chip_seq %>% filter(V9=="SRR577356_SRR577357.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me3_2_name = V1, H3K4me3_2_start = V2, H3K4me3_2_end = V3, H3K4me3_2_RPKM = V10) %>% select(H3K4me3_2_RPKM)

h3k36me3  = chip_seq %>% filter(V9=="SRR577400_SRR577401.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K36me3_name = V1, H3K36me3_start = V2, H3K36me3_end = V3, H3K36me3_RPKM = V10) %>% select(H3K36me3_RPKM)

h3k36me3_2  = chip_seq %>% filter(V9=="SRR227435_SRR227436.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K36me3_2_name = V1, H3K36me3_2_start = V2, H3K36me3_2_end = V3, H3K36me3_2_RPKM = V10) %>% select(H3K36me3_2_RPKM)

h3k27me3  = chip_seq %>% filter(V9=="SRR577370_SRR577371.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K27me3_name = V1, H3K27me3_start = V2, H3K27me3_end = V3, H3K27me3_RPKM = V10) %>% select(H3K27me3_RPKM)

h3k27me3_2  = chip_seq %>% filter(V9=="SRR350910_SRR227607_SRR227608.fastq") %>% select(V1,V2,V3,V10) %>% mutate(H3K27me3_2_name = V1, H3K27me3_2_start = V2, H3K27me3_2_end = V3, H3K27me3_2_RPKM = V10) %>% select(H3K27me3_2_RPKM)

h3k4me1  = chip_seq %>% filter(V8=="H3K4me1") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me1_name = V1, H3K4me1_start = V2, H3K4me1_end = V3, H3K4me1_RPKM = V10) %>% select(H3K4me1_RPKM)

h2afz  = chip_seq %>% filter(V8=="H2AFZ") %>% select(V1,V2,V3,V10) %>% mutate(H2AFZ_name = V1, H2AFZ_start = V2, H2AFZ_end = V3, H2AFZ_RPKM = V10) %>% select(H2AFZ_RPKM)

h3k27ac  = chip_seq %>% filter(V8=="H3K27ac") %>% select(V1,V2,V3,V10) %>% mutate(H3K27ac_name = V1, H3K27ac_start = V2, H3K27ac_end = V3, H3K27ac_RPKM = V10) %>% select(H3K27ac_RPKM)

h3k4me2  = chip_seq %>% filter(V8=="H3K4me2") %>% select(V1,V2,V3,V10) %>% mutate(H3K4me2_name = V1, H3K4me2_start = V2, H3K4me2_end = V3, H3K4me2_RPKM = V10) %>% select(H3K4me2_RPKM)

h3k79me2  = chip_seq %>% filter(V8=="H3K79me2") %>% select(V1,V2,V3,V10) %>% mutate(H3K79me2_name = V1, H3K79me2_start = V2, H3K79me2_end = V3, H3K79me2_RPKM = V10) %>% select(H3K79me2_RPKM)

h3k9ac  = chip_seq %>% filter(V8=="H3K9ac") %>% select(V1,V2,V3,V10) %>% mutate(H3K9ac_name = V1, H3K9ac_start = V2, H3K9ac_end = V3, H3K9ac_RPKM = V10) %>% select(H3K9ac_RPKM)

h3k9me3 = chip_seq %>% filter(V8=="H3K9me3") %>% select(V1,V2,V3,V10) %>% mutate(H3K9me3_name = V1, H3K9me3_start = V2, H3K9me3_end = V3, H3K9me3_RPKM = V10) %>% select(H3K9me3_RPKM)

h4k20 = chip_seq %>% filter(V8=="H4K20me1") %>% select(V1,V2,V3,V10) %>% mutate(H4K20me1_name = V1, H4K20me1_start = V2, H4K20me1_end = V3, H4K20me1_RPKM = V10) %>% select(H4K20me1_RPKM)

############################################################################################################################
# Filtered DNAse I hypersenstivity sites for HeLa cells:

DNAse_I = DNaseI %>% filter(V8 == "DNAse-seq") %>% select(V1,V2,V3,V10) %>% mutate(DNAse_name = V1, DNAse_start = V2, DNAse_end = V3, DNAse_RPKM = V10) %>% select(DNAse_RPKM)

############################################################################################################################
# Filter XR-seq Data for GM12878 cells:

XR_1 = xr_seq %>% filter(V9 == "Cisplatin_Replicate1") %>% select(V1, V2, V3, V10) %>% mutate(xr1_name = V1, xr1_start = V2, xr1_end = V3, xr1_RPKM = V10) %>% select(xr1_name, xr1_start, xr1_end, xr1_RPKM)

XR_2 = xr_seq %>% filter(V9 == "Cisplatin_Replicate2") %>% select(V1, V2, V3, V10) %>% mutate(xr2_name = V1, xr2_start = V2, xr2_end = V3, xr2_RPKM = V10) %>% select(xr2_RPKM)
############################################################################################################################
# Filter XR-seq Data for HeLa cells:

XR_1 = xr_seq %>% filter(V9 == "HXA64A1_ATCACG") %>% select(V1, V2, V3, V10) %>% mutate(xr1_name = V1, xr1_start = V2, xr1_end = V3, xr1_RPKM = V10) %>% select(xr1_name, xr1_start, xr1_end, xr1_RPKM)

XR_2 = xr_seq %>% filter(V9 == "HXA64B7_CAGATC") %>% select(V1, V2, V3, V10) %>% mutate(xr2_name = V1, xr2_start = V2, xr2_end = V3, xr2_RPKM = V10) %>% select(xr2_RPKM)

XR_3 = xr_seq %>% filter(V9 == "HXACA4_TGACCA") %>% select(V1, V2, V3, V10) %>% mutate(xr3_name = V1, xr3_start = V2, xr3_end = V3, xr3_RPKM = V10) %>% select(xr3_RPKM)

XR_4 = xr_seq %>% filter(V9 == "HXACB10_TAGCTT") %>% select(V1, V2, V3, V10) %>% mutate(xr4_name = V1, xr4_start = V2, xr4_end = V3, xr4_RPKM = V10) %>% select(xr4_RPKM)
############################################################################################################################
# Filter Damage-seq Data for GM12878 cells:

dmg_1 = dmg_seq %>% filter(V9 == "Cisplatin_Replicate1") %>% select(V1, V2, V3, V10) %>% mutate(dmg1_name = V1, dmg1_start = V2, dmg1_end = V3, dmg1_RPKM = V10) %>% select(dmg1_RPKM)

dmg_2 = dmg_seq %>% filter(V9 == "Cisplatin_Replicate2") %>% select(V1, V2, V3, V10) %>% mutate(dmg2_name = V1, dmg2_start = V2, dmg2_end = V3, dmg2_RPKM = V10) %>% select(dmg2_RPKM)
############################################################################################################################
# Filter Damage-seq Data for HeLa cells:

dmg_1 = dmg_seq %>% filter(V9 == "HDA64A1_ATCACG") %>% select(V1, V2, V3, V10) %>% mutate(dmg1_name = V1, dmg1_start = V2, dmg1_end = V3, dmg1_RPKM = V10) %>% select(dmg1_RPKM)

dmg_2 = dmg_seq %>% filter(V9 == "HDA64B19_GTGAAA") %>% select(V1, V2, V3, V10) %>% mutate(dmg2_name = V1, dmg2_start = V2, dmg2_end = V3, dmg2_RPKM = V10) %>% select(dmg2_RPKM)

dmg_3 = dmg_seq %>% filter(V9 == "HDACA6_GCCAAT") %>% select(V1, V2, V3, V10) %>% mutate(dmg3_name = V1, dmg3_start = V2, dmg3_end = V3, dmg3_RPKM = V10) %>% select(dmg3_RPKM)

dmg_4 = dmg_seq %>% filter(V9 == "HDACB23_GAGTGG") %>% select(V1, V2, V3, V10) %>% mutate(dmg4_name = V1, dmg4_start = V2, dmg4_end = V3, dmg4_RPKM = V10) %>% select(dmg4_RPKM)

############################################################################################################################
rm(xr_seq,dmg_seq,chip_seq,DNaseI)
# Create Data Frame:

# DF for HeLa Cells:
#df_allseq = data.frame(XR_1, XR_2, XR_3, XR_4, dmg_1, dmg_2, dmg_3, dmg_4, h4k20, h3k9me3, h3k9ac, h3k36me3, h3k79me2, h2afz, h3k4me2, h3k27ac, h3k27me3, h3k4me1, h3k4me3, h3k4me3_2, h3k4me3_3,h3k27me3_2, h3k36me3_2, DNAse_I)
# DF for GM12878 Cells:
df_allseq = data.frame(XR_1, XR_2, dmg_1, dmg_2, h4k20, h3k9me3, h3k9ac, h3k36me3, h3k79me2, h2afz, h3k4me2, h3k27ac, h3k27me3, h3k4me1, h3k4me3, h3k4me3_2, h3k27me3_2, h3k36me3_2)

head(df_allseq)

############################################################################################################################
# Write down dataframe to a file:
setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Stock_Seq_Files/Only_RPKM/Subsampled')

write.table(df_allseq, file="STOCK_allseq_subsampled_RPKM.bed", quote=F, sep="\t", row.names=F, col.names=TRUE)
