library(tidyverse)
library(magrittr)
library(GGally)

##  Open the big three

setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Stock_Seq_Files')

repair = read.table('total_xr-seqs.bed',header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(repair)

chipseq = read.table('total_chip-seq.bed',header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)
head(chipseq)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Damage-seq_Data/No_ControlandSignDiversity_Data')
damage =  read.table('dmg_seqs_no_signs-woHeader.csv',header = FALSE, sep="\t", stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(damage)    ## Strands needed to be combined by adding corresponding RPKMs of bot  plus and minus rows. 

## Filtereing and Selection Part:

chr_damage_seq = damage %>% filter(V9 == "Cisplatin_Replicate1") %>% select(V1, V2, V3) %>% mutate(chr_name = V1, chr_start = V2, chr_end = V3) %>% select(chr_name, chr_start, chr_end)
head(chr_damage_seq)

damage_cis1 = damage %>% filter(V9 == "Cisplatin_Replicate1") %>% select(V10) %>% mutate(dmg_1 = V10) %>% select(dmg_1)
damage_cis2 = damage %>% filter(V9 == "Cisplatin_Replicate2") %>% select(V10) %>% mutate(dmg_2 = V10) %>% select(dmg_2)
head(damage_cis1)
head(damage_cis2)

repair_cis1 = repair %>% filter(V9 == "Cisplatin_Replicate1") %>% select(V10) %>% mutate(xr_1 = V10) %>% select(xr_1)
repair_cis2 = repair %>% filter(V9 == "Cisplatin_Replicate2") %>% select(V10) %>% mutate(xr_2 = V10) %>% select(xr_2)
head(repair_cis1)
head(repair_cis2)

h4k20 = chipseq %>% filter(V8 == "H4K20me1") %>% select(V10) %>% mutate(RPKM_H4K20me1 = V10) %>% select(RPKM_H4K20me1)
h3k9me3 = chipseq %>% filter(V8 == "H3K9me3") %>% select(V10) %>% mutate(RPKM_H3K9me3 = V10) %>% select(RPKM_H3K9me3)
h3k9ac = chipseq %>% filter(V8 == "H3K9ac") %>% select(V10) %>% mutate(RPKM_H3K9ac = V10) %>% select(RPKM_H3K9ac)
h3k36_1 = chipseq %>% filter(V9 == "SRR227435_SRR227436.fastq") %>% select(V10) %>% mutate(RPKM_H3K36me3_1 = V10) %>% select(RPKM_H3K36me3_1)
h3k36_2 = chipseq %>% filter(V9 == "SRR577400_SRR577401.fastq") %>% select(V10) %>% mutate(RPKM_H3K36me3_2 = V10) %>% select(RPKM_H3K36me3_2)
h3k79me2 = chipseq %>% filter(V8 == "H3K79me2") %>% select(V10) %>% mutate(RPKM_H3K79me2 = V10) %>% select(RPKM_H3K79me2)
h2afz = chipseq %>% filter(V8 == "H2AFZ") %>% select(V10) %>% mutate(RPKM_H2AFZ = V10) %>% select(RPKM_H2AFZ)
h3k4me2 = chipseq %>% filter(V8 == "H3K4me2") %>% select(V10) %>% mutate(RPKM_H3K4me2 = V10) %>% select(RPKM_H3K4me2)
h3k27ac = chipseq %>% filter(V8 == "H3K27ac") %>% select(V10) %>% mutate(RPKM_H3K27ac = V10) %>% select(RPKM_H3K27ac)
h3k27_1 = chipseq %>% filter(V9 == "SRR350910_SRR227607_SRR227608.fastq") %>% select(V10) %>% mutate(RPKM_H3K27me3_1 = V10) %>% select(RPKM_H3K27me3_1)
h3k27_2 = chipseq %>% filter(V9 == "SRR577370_SRR577371.fastq") %>% select(V10) %>% mutate(RPKM_H3K27me3_2 = V10) %>% select(RPKM_H3K27me3_2)
h3k4me1 = chipseq %>% filter(V8 == "H3K4me1") %>% select(V10) %>% mutate(RPKM_H3K4me1 = V10) %>% select(RPKM_H3K4me1)
h3k4me3_1 = chipseq %>% filter(V9 == "SRR5331171_SRR5331172.fastq") %>% select(V10) %>% mutate(RPKM_H3K4me3_1 = V10) %>% select(RPKM_H3K4me3_1)
h3k4me3_2 = chipseq %>% filter(V9 == "SRR577356_SRR577357.fastq") %>% select(V10) %>% mutate(RPKM_H3K4me3_2 = V10) %>% select(RPKM_H3K4me3_2)

head(h3k36_1)
head(h3k36_2)
head(h3k27_1)
head(h3k27_2)
head(h3k4me3_1)
head(h3k4me3_2)

rm(repair, damage, chipseq)

# Data Frame Creator:
df_allseq <- data.frame(chr_damage_seq,repair_cis1, repair_cis2, damage_cis1, damage_cis2, h4k20, h3k9me3, h3k9ac, h3k36_1, h3k79me2, h2afz, h3k4me2, h3k27ac, h3k27_1, h3k4me1, h3k4me3_1, h3k4me3_2, h3k27_2, h3k36_2)
head(df_allseq)

##  For saving file as bed
setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/GM12878-Cells/Stock_Seq_Files/Only_RPKM')
write.table(df_allseq, "STOCK_allseq_RPKM.bed",sep="\t",row.names=FALSE)
