#   Bu kod her türlü seq file'ın "0" RPKM değerine sahip satırlarını atmaya ve repairi damage'a bölmeye yarar.
#   #   Sütun numaraları farklılık gösterebilir.
#   İlk 3 sütun hariç geriye kalan tüm sütunlar RPKM değeri olmalı!
#   Remove rows with respect to corresponding Damage-seq
"""
stock_file_path = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//xr_damage_chip_DNase-seq_Data//Stock_Only_RPKM_Data//'
stock_file = open(stock_file_path + 'STOCK_allseq_RPKM_woHEADER.bed','r')

##  Row removal and Normalization for 6-4PP Damage ##
without_zero_rows_PPD_RepA = []
without_zero_rows_PPD_RepB = []
for indiv_lines in stock_file:
    splitted_lines = indiv_lines.split('\t')
    PPD_three_columns = splitted_lines[0:3]         #   CHR name, start and end coordinates
    PPD_histone_markers = splitted_lines[11:26]     #   The RPKM values column of ChIP-seq histone markers.
    XR_seq_PPD_RepA = float(splitted_lines[3])      #   The RPKM values column of XR-seq Replicate A.
    XR_seq_PPD_RepB = float(splitted_lines[4])      #   The RPKM values column of XR-seq Replicate B.
    DMG_PPD_column_RepA = float(splitted_lines[7])  #   The Damage-seq column of Replicate A of 6-4PP RepA damage.
    DMG_PPD_column_RepB = float(splitted_lines[8])  #   The Damage-seq column of Replicate B of 6-4PP RepB damage.
    PPD_DNase_column = float(splitted_lines[26])
    if DMG_PPD_column_RepA != 0:
        R1_one = float(XR_seq_PPD_RepA / DMG_PPD_column_RepA)       #   This part directly adds rows NOT having 0 RPKM into list so called without_zero_rows_PPD_RepA.
        without_zero_rows_PPD_RepA.append([str(PPD_three_columns),float(R1_one),str(PPD_histone_markers),str(PPD_DNase_column)])
    if DMG_PPD_column_RepB != 0:
        R1_two = float(XR_seq_PPD_RepB / DMG_PPD_column_RepB)       # This part directly adds rows NOT having 0 RPKM into list so called without_zero_rows_PPD_RepB.  
        without_zero_rows_PPD_RepB.append([str(PPD_three_columns),float(R1_two),str(PPD_histone_markers),str(PPD_DNase_column)])

print("Number of elements in without_zero_rows_PPD_RepA list: " , len(without_zero_rows_PPD_RepA))
print("Number of elements in without_zero_rows_PPD_RepB list: " , len(without_zero_rows_PPD_RepB))
stock_file.close()
"""
stock_file_path = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//xr_damage_chip_DNase-seq_Data//Stock_Only_RPKM_Data//'
stock_file = open(stock_file_path + 'STOCK_allseq_RPKM_woHEADER.bed','r')

##  Row removal and Normalization for 6-4PP Damage ##
without_zero_rows_CPD_RepA = []
without_zero_rows_CPD_RepB = []
for indiv_lines_2 in stock_file:
    splitted_lines_2 = indiv_lines_2.split('\t')
    CPD_three_columns = splitted_lines_2[0:3]         #   CHR name, start and end coordinates
    CPD_histone_markers = splitted_lines_2[11:26]     #   The RPKM values column of ChIP-seq histone markers.
    XR_seq_CPD_RepA = float(splitted_lines_2[5])      #   The RPKM values column of XR-seq Replicate A.
    XR_seq_CPD_RepB = float(splitted_lines_2[6])      #   The RPKM values column of XR-seq Replicate B.
    DMG_CPD_column_RepA = float(splitted_lines_2[9])  #   The Damage-seq column of Replicate A of 6-4PP RepA damage.
    DMG_CPD_column_RepB = float(splitted_lines_2[10])  #   The Damage-seq column of Replicate B of 6-4PP RepB damage.
    CPD_DNase_column = float(splitted_lines_2[26])
    if DMG_CPD_column_RepA != 0:
        R2_one = float(XR_seq_CPD_RepA / DMG_CPD_column_RepA)       #   This part directly adds rows NOT having 0 RPKM into list so called without_zero_rows_CPD_RepA.
        without_zero_rows_CPD_RepA.append([str(CPD_three_columns),float(R2_one),str(CPD_histone_markers),str(CPD_DNase_column)])
    if DMG_CPD_column_RepB != 0:
        R2_two = float(XR_seq_CPD_RepB / DMG_CPD_column_RepB)       # This part directly adds rows NOT having 0 RPKM into list so called without_zero_rows_CPD_RepB.  
        without_zero_rows_CPD_RepB.append([str(CPD_three_columns),float(R2_two),str(CPD_histone_markers),str(CPD_DNase_column)])

print("Number of elements in without_zero_rows_CPD_RepA list: " , len(without_zero_rows_CPD_RepA))
print("Number of elements in without_zero_rows_CPD_RepB list: " , len(without_zero_rows_CPD_RepB))
stock_file.close()

##################################################### CREATE FILE #####################################################

######  Write down the remaining of zero removed rows into a file.  #####
#HDA64A1_ATCACG                     HDA64B19_GTGAAA                     HDACA6_GCCAAT                       HDACB23_GAGTGG
#without_zero_rows_PPD_RepA     without_zero_rows_PPD_RepB      without_zero_rows_CPD_RepA              without_zero_rows_CPD_RepB

total = 0
directory = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//xr_damage_chip_DNase-seq_Data//Normal Filtered_Normalized_allseq-Data//'
with open(directory + 'Zero_removed_divided_HDACB23_GAGTGG_normalized_allseqs_w_ChrStartEnd.csv','w') as df:
    for wr in without_zero_rows_CPD_RepB:
        df.write('{}\n'.format(wr))       #       Do not put \n front of '{}'.
        total += 1
print("The newly generated file is having " , total , ' number of rows!')
########################################################################################################################
#   The newly created Datium Structure will be:   #
#   Chr names   Chr Start   Chr END     6-4PP (XR-seq RepA / DMG RepA)    ChIP-seq (15 columns)
#   Chr names   Chr Start   Chr END     6-4PP(XR-seq RepB / DMG RepB)    ChIP-seq (15 columns)
#   Chr names   Chr Start   Chr END     CPD(XR-seq RepA / DMG RepA)    ChIP-seq (15 columns)
#   Chr names   Chr Start   Chr END     CPD(XR-seq RepB / DMG RepB)    ChIP-seq (15 columns)
