library(tidyverse)
library(magrittr)
library(GGally)

setwd('C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells/ReplicateA-Normalized-woZero-Divided-allseq-Data')

PP_repA_seqs <- read.csv("Zero_removed_divided_HDA64A1_ATCACG_normalized_allseqs_w_ChrStartEnd.xlsx", header = FALSE, 
                      sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

CPD_repA_seqs <- read.csv("Zero_removed_divided_HDACA6_GCCAAT_normalized_allseqs_w_ChrStartEnd.xlsx", header = FALSE, 
                           sep=",",stringsAsFactors=FALSE, quote="\"", na.strings = "NA", fill=TRUE)

head(PP_repA_seqs)
head(CPD_repA_seqs)

PP_repA_seqs_filtered = PP_repA_seqs %>% select(V1:V20) %>% mutate(Start_PP = as.factor(V2), End_PP = as.factor(V3),
                                                                                   XR_64PP_seq_RepA = as.numeric(V4),
                                                                                   XR_64PP_seq_RepB = as.numeric(V5),
                                                                                   Histones_64PP = as.numeric(V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20),
                                                                                   CHR_PP = as.factor(V1)) %>% select(CHR_PP, Start_PP, End_PP, XR_64PP_seq_RepA, XR_64PP_seq_RepB, Histones_64PP)

CPD_repA_seqs_filtered = CPD_repA_seqs %>% select(V1:V20) %>% mutate(Start_CPD = as.factor(V2), End_CPD = as.factor(V3),
                                                                     XR_CPD_seq_RepA = as.numeric(V4),
                                                                     XR_CPD_seq_RepB = as.numeric(V5),
                                                                     Histones_CPD = as.numeric(V6,V7,V8,V9,V10,V11,V12,V13,V14,V15,V16,V17,V18,V19,V20),
                                                                     CHR_CPD = as.factor(V1)) %>% select(CHR_CPD, Start_CPD, End_CPD, XR_CPD_seq_RepA, XR_CPD_seq_RepB, Histones_CPD)

HDA64A1_ATCACGc = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% mutate(HDA64A1_ATCACG = XR_64PP_seq_RepA) %>% 
  select(HDA64A1_ATCACG, CHR_PP, Start_PP, End_PP, Histones_64PP)

HDA64B19_GTGAAAc = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% mutate(HDA64B19_GTGAAA = XR_64PP_seq_RepB) %>% 
  select(HDA64B19_GTGAAA, CHR_PP, Start_PP, End_PP, Histones_64PP)

HDACA6_GCCAATc = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% mutate(HDACA6_GCCAAT = XR_CPD_seq_RepA) %>% 
  select(HDACA6_GCCAAT, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)

HDACB23_GAGTGGc = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% mutate(HDACB23_GAGTGG = XR_CPD_seq_RepB) %>% 
  select(HDACB23_GAGTGG, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)

# Data frame
HDA64A1_ATCACGc_sorted = HDA64A1_ATCACGc %>% arrange(desc(HDA64A1_ATCACG))

HDA64B19_GTGAAAc_sorted = HDA64B19_GTGAAAc %>% arrange(desc(HDA64B19_GTGAAA))

HDACA6_GCCAATc_sorted = HDACA6_GCCAATc %>% arrange(desc(HDACA6_GCCAAT))

HDACB23_GAGTGGc_sorted = HDACB23_GAGTGGc %>% arrange(desc(HDACB23_GAGTGG))

head(HDA64A1_ATCACGc_sorted)
#tail(HDA64A1_ATCACGc_sorted)

head(HDA64B19_GTGAAAc_sorted)
#tail(HDA64B19_GTGAAAc_sorted)

head(HDACA6_GCCAATc_sorted)
#tail(HDACA6_GCCAATc_sorted)

head(HDACB23_GAGTGGc_sorted)
#tail(HDACB23_GAGTGGc_sorted)
######################################################

HDA64A1_ATCACGc_tail = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepA == "1") %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP)
HDA64A1_ATCACGc_tail_2 = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepA == "2") %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP)
HDA64A1_ATCACGc_tail_3 = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepA == "3") %>% select(XR_64PP_seq_RepA, CHR_PP, Start_PP, End_PP, Histones_64PP)


HDA64B19_GTGAAAc_tail = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepB == "1") %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP)
HDA64B19_GTGAAAc_tail_2 = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepB == "2") %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP)
HDA64B19_GTGAAAc_tail_3 = PP_repA_seqs_filtered %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP) %>% filter(XR_64PP_seq_RepB == "3") %>% select(XR_64PP_seq_RepB, CHR_PP, Start_PP, End_PP, Histones_64PP)

HDACA6_GCCAATc_tail = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepA == "1") %>%  select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)
HDACA6_GCCAATc_tail_2 = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepA == "2") %>%  select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)
HDACA6_GCCAATc_tail_3 = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepA == "3") %>%  select(XR_CPD_seq_RepA, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)


HDACB23_GAGTGGc_tail = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepB == "1") %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)
HDACB23_GAGTGGc_tail_2 = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepB == "2") %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)
HDACB23_GAGTGGc_tail_3 = CPD_repA_seqs_filtered %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD) %>% filter(XR_CPD_seq_RepB == "3") %>% select(XR_CPD_seq_RepB, CHR_CPD, Start_CPD, End_CPD, Histones_CPD)

# Data frame for smallest repair
HDA64A1_ATCACGc_tail_sort = HDA64A1_ATCACGc_tail %>% arrange(XR_64PP_seq_RepA)
HDA64A1_ATCACGc_tail_sort_2 = HDA64A1_ATCACGc_tail_2 %>% arrange(XR_64PP_seq_RepA)
HDA64A1_ATCACGc_tail_sort_3 = HDA64A1_ATCACGc_tail_3 %>% arrange(XR_64PP_seq_RepA)

HDA64B19_GTGAAAc_tail_sort = HDA64B19_GTGAAAc_tail %>% arrange(XR_64PP_seq_RepB)
HDA64B19_GTGAAAc_tail_sort_2 = HDA64B19_GTGAAAc_tail_2 %>% arrange(XR_64PP_seq_RepB)
HDA64B19_GTGAAAc_tail_sort_3 = HDA64B19_GTGAAAc_tail_3 %>% arrange(XR_64PP_seq_RepB)

HDACA6_GCCAATc_tail_sort = HDACA6_GCCAATc_tail %>% arrange(XR_CPD_seq_RepA)
HDACA6_GCCAATc_tail_sort_2 = HDACA6_GCCAATc_tail_2 %>% arrange(XR_CPD_seq_RepA)
HDACA6_GCCAATc_tail_sort_3 = HDACA6_GCCAATc_tail_3 %>% arrange(XR_CPD_seq_RepA)

HDACB23_GAGTGGc_tail_sort = HDACB23_GAGTGGc_tail %>% arrange(XR_CPD_seq_RepB)
HDACB23_GAGTGGc_tail_sort_2 = HDACB23_GAGTGGc_tail_2 %>% arrange(XR_CPD_seq_RepB)
HDACB23_GAGTGGc_tail_sort_3 = HDACB23_GAGTGGc_tail_3 %>% arrange(XR_CPD_seq_RepB)

head(HDA64A1_ATCACGc_tail_sort)
head(HDA64A1_ATCACGc_tail_sort_2)
head(HDA64A1_ATCACGc_tail_sort_3)

head(HDA64B19_GTGAAAc_tail_sort)
head(HDA64B19_GTGAAAc_tail_sort_2)
head(HDA64B19_GTGAAAc_tail_sort_3)

head(HDACA6_GCCAATc_tail_sort)
head(HDACA6_GCCAATc_tail_sort_2)
head(HDACA6_GCCAATc_tail_sort_3)

head(HDACB23_GAGTGGc_tail_sort)
head(HDACB23_GAGTGGc_tail_sort_2)
head(HDACB23_GAGTGGc_tail_sort_3)
##  Create File ##
setwd("C:/Users/Arda/Documents/TEZ/TEZ-Data/HeLa-Cells")

##  New file will have HEADER so there will be 606224 rows
write.table(data_frame_generator, "XR_over_DMG_all-seq-including_Peak_Data.csv", sep = "\t", row.names = FALSE)