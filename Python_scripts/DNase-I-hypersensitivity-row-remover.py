##          This script removes rows from DNase I hypersensitivity regions with respect to Damage-seq file.
path_DMG = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//Peak-Repair-DMG-ChIP-seq_wo_Cntrl-DATA//Only_RPKM_Data//STOCK-Data_wo_Cntrl//WithCHRStartEnd//'

DMG_file = open(path_DMG + 'XR-DMG-ChIP_wo_Cntrl-Peak_Data_w_ChrStartEnd.csv','r')

DMG_1_list = []
DMG_2_list = []
DMG_3_list = []
DMG_4_list = []
for DMG_line in DMG_file:                       #   This part select columns from Damage-seq files having RPKM and putting them into list.
    splitted = DMG_line.split("\t")             #   For file C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//Peak-Repair-DMG-ChIP-seq_wo_Cntrl-DATA//Only_RPKM_Data//STOCK-Data_wo_Cntrl//WithCHRStartEnd//'XR-DMG-ChIP_wo_Cntrl-Peak_Data_w_ChrStartEnd.csv
    PPD_RepA_column = splitted[7]               #   #   PPD_RepA -> column 7
    DMG_1_list.append(PPD_RepA_column)
    PPD_RepB_column = splitted[8]               #   #   PPD_RepB -> column 8
    DMG_2_list.append(PPD_RepB_column)
    CPD_RepA_column = splitted[9]               #   #   CPD_RepA -> column 9
    DMG_3_list.append(CPD_RepA_column)
    CPD_RepB_column = splitted[10]              #   #   CPD_RepB -> column 10
    DMG_4_list.append(CPD_RepB_column)

print("The number of RPKMs within PPD_RepA list: " , len(DMG_1_list))
print("The number of RPKMs within PPD_RepB list: " , len(DMG_2_list))
print("The number of RPKMs within CPD_RepA list: " , len(DMG_3_list))
print("The number of RPKMs within CPD_RepB list: " , len(DMG_4_list))
DMG_file.close()

#############################################################################################################################################3
path_DNase = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//DNaseI_Hypersensitivity_Sites//'

DNase_file = open(path_DNase + '99_SRR352412_SRR352413_SRR352414.fastq_RPKM_added.bed','r')

DNase_I_hypersensitivity_list = []          #   This part selects whole row from DNase I hypersensitivity regions file and put them one by one into list.
for DNase_line in DNase_file:
    DNase_I_hypersensitivity_list.append(DNase_line)

print('The number of RPKM elements found in DNase_I_hypersensitivity: ' , len(DNase_I_hypersensitivity_list))
DNase_file.close()

#   Here using the Damage seq file, DNase I hypersensitivity rows are removed and remaining rows are put into another list together with all BED format information.
#   #   There are 4 lists due to presence of four different damage files therefore, there will be 4 different output DNase-seq files.
zero_removed_DNaseI_PPD_RepA = []
zero_removed_DNaseI_PPD_RepB = []
zero_removed_DNaseI_CPD_RepA = []
zero_removed_DNaseI_CPD_RepB = []
for k in range(0,len(DNase_I_hypersensitivity_list)):
    if DMG_1_list[k] != "0":
        zero_removed_DNaseI_PPD_RepA.append(DNase_I_hypersensitivity_list[k])
    if DMG_2_list[k] != "0":
        zero_removed_DNaseI_PPD_RepB.append(DNase_I_hypersensitivity_list[k])
    if DMG_3_list[k] != "0":
        zero_removed_DNaseI_CPD_RepA.append(DNase_I_hypersensitivity_list[k])
    if DMG_4_list[k] != "0":
        zero_removed_DNaseI_CPD_RepB.append(DNase_I_hypersensitivity_list[k])

print('The number of elements within list zero_removed_DNaseI_PPD_RepA: ' , len(zero_removed_DNaseI_PPD_RepA))
print('The number of elements within list zero_removed_DNaseI_PPD_RepB: ' , len(zero_removed_DNaseI_PPD_RepB))
print('The number of elements within list zero_removed_DNaseI_CPD_RepA: ' , len(zero_removed_DNaseI_CPD_RepA))
print('The number of elements within list zero_removed_DNaseI_CPD_RepB: ' , len(zero_removed_DNaseI_CPD_RepB))
"""
##  #   Write file #    ##

#   HDA64A1_ATCACG                      HDA64B19_GTGAAA                         HDACA6_GCCAAT                       HDACB23_GAGTGG
#   zero_removed_DNaseI_PPD_RepA      zero_removed_DNaseI_PPD_RepB       zero_removed_DNaseI_CPD_RepA           zero_removed_DNaseI_CPD_RepB
total = 0
directory = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//DNaseI_Hypersensitivity_Sites//Zero_rows//'
with open(directory + 'Zero_removed_HDA64A1_ATCACG_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv','w') as df:
    for wr in zero_removed_DNaseI_PPD_RepA:
        df.write('{}\n'.format(wr))       #       Do not put \n front of '{}'.
        total += 1
print("The newly generated file is having " , total , ' number of rows!')
"""

### REMOVE ROWS WITH RESPECT TO DAMAGE-seq REPLICATE A  ###
###############################################################################################################
stock_file_path = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//Peak-Repair-DMG-ChIP-seq_wo_Cntrl-DATA//Only_RPKM_Data//STOCK-Data_wo_Cntrl//WithCHRStartEnd//'
stock_file = open(stock_file_path + 'XR-DMG-ChIP_wo_Cntrl-Peak_Data_w_ChrStartEnd.csv','r')

##  Row removal and Normalization of DNase I hypersensitivity region file from 6-4PP Damage ##
PPD_DMG_RepA_list = []
for indiv_lines in stock_file:
    splitted_lines = indiv_lines.split('\t')
    DMG_PPD_column = float(splitted_lines[7])
    PPD_DMG_RepA_list.append(DMG_PPD_column)        #   Append every 6-4 DMG replicate A rows into a list (PPD_DMG_RepA_list).
print("Number of elements within list PPD_DMG_RepA: " , len(PPD_DMG_RepA_list))

stock_file.close()

###############################################################################################################
stock_file = open(stock_file_path + 'XR-DMG-ChIP_wo_Cntrl-Peak_Data_w_ChrStartEnd.csv','r')

##  Row removal and Normalization of DNase I hypersensitivity region from file CPD Damage    ##
CPD_DMG_RepA_list = []
for indiv_lines_2 in stock_file:
    splitted_lines_2 = indiv_lines_2.split('\t')
    DMG_CPD_column = float(splitted_lines_2[9])
    CPD_DMG_RepA_list.append(DMG_CPD_column)        #   Append every CPD DMG replicate A rows into a list (CPD_DMG_RepA_list)
print("Number of elements within list CPD_DMG_RepA: " , len(CPD_DMG_RepA_list))

stock_file.close()

###############################################################################################################
##  This part selects whole row from DNase I hypersensitivity regions file and put them one by one into list.
path_DNase = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//DNaseI_Hypersensitivity_Sites//'

DNase_file = open(path_DNase + '99_SRR352412_SRR352413_SRR352414.fastq_RPKM_added.bed','r')

DNase_I_hypersensitivity_list = []
for DNase_line in DNase_file:
    DNase_I_hypersensitivity_list.append(DNase_line)

print('The number of RPKM elements found in DNase_I_hypersensitivity: ' , len(DNase_I_hypersensitivity_list))
DNase_file.close()

##   Remove rows of DNase I hypersensitivity with respect to 6-4PP damage replicate A zero rows.
rows_removed_wr_PPD_RepA_list = []
for n in range(0,len(DNase_I_hypersensitivity_list)):
    if PPD_DMG_RepA_list[n] != 0:
        rows_removed_wr_PPD_RepA_list.append(DNase_I_hypersensitivity_list[n])
print('The number of elements within list rows_removed_wr_PPD_RepA_list: ' , len(rows_removed_wr_PPD_RepA_list))

##   Remove rows of DNase I hypersensitivity with respect to CPD damage replicate A zero rows.
rows_removed_wr_CPD_RepA_list = []
for k in range(0,len(DNase_I_hypersensitivity_list)):
    if CPD_DMG_RepA_list[k] != 0:
        rows_removed_wr_CPD_RepA_list.append(DNase_I_hypersensitivity_list[k])
print('The number of elements within list rows_removed_wr_CPD_RepA_list: ' , len(rows_removed_wr_CPD_RepA_list))

#####################################################################################################################
"""
##  #   Write file #    ##

#   HDA64A1_ATCACG   HDA64B19_GTGAAA        HDACA6_GCCAAT  HDACB23_GAGTGG
#       rows_removed_wr_PPD_RepA_list     rows_removed_wr_CPD_RepA_list           
total = 0
directory = 'C://Users//Arda//Documents//TEZ//TEZ-Data//HeLa-Cells//DNaseI_Hypersensitivity_Sites//Zero_rows//ReplicateA_Normalized//'
with open(directory + 'Zero_removed_HDACA6_GCCAAT_normalized_DNse_I_hypersensitivity_w_ChrStartEnd.csv','w') as df:
    for wr in rows_removed_wr_CPD_RepA_list:
        df.write('{}\n'.format(wr))       #       Do not put \n front of '{}'.
        total += 1
print("The newly generated file is having " , total , ' number of rows!')
"""
