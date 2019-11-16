    ##### This script is written for selecting desired column from txt or BED file. #####
    ### Script gives a list containing list of elements of desired column's each row. ###

"CAUTION! --Files should be equally tab or space-delimited--"
import numpy as np
import io

column_elements_list = [""]
with open('BED_2.txt','r') as df:
    for loop in df:
        split_tab = loop.split("\t")
        #print(split_tab)                   #   split_tab creates a list.
        filtering_column  = split_tab[9]
        #print(filtering_column)
        #print(type(filtering_column))      #   filtering_column is a string not list.
        desired_column = filtering_column.replace("\n","")  #   The new line escape operator "\n" is removed.
        column_elements_list.append(desired_column)
    print(column_elements_list)             # Now! The list so called "column_elements_list" has the desired column of elements.