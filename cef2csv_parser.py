# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 13:39:40 2018

@author: Choon
"""


import os
import io

cef_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
#cef_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"


inputfile = open(os.getcwd()+"\\"+cef_file_name+".cef")
text = inputfile.readlines()


i=0
while i<len(text):
    if ((text[i]).split())[0] == 'FILE_NAME':
        if ((text[i]).split())[2][-5:] == '.cef"':
            file_name = ((text[i]).split())[2][1:-5]
            break
        else:
            raise Exception("Check .cef file format")
    i+=1


var_arr = []
while i<len(text):
    if len(((text[i]).split())) > 0:
        if ((text[i]).split())[0] == 'START_VARIABLE':
            var_arr.append(((text[i]).split())[2])
    if text[i] == 'DATA_UNTIL = EOF\n':
        break
    i+=1
    
data_start_pos = i+1

#########################################
csv_file = io.open(os.getcwd() +"\\"+ file_name + ".csv",'w')

data_arr=[]

i=data_start_pos
while ((text[i]).split())[0] != '!RECORDS=':
    data_row_str = ((text[i]).split())[0]
    csv_file.write(data_row_str+"\n")
    
    i+=1








    

#ceflib.read(os.getcwd()+"\\C1_CP_FGM_SPIN__20180110_000000_20180117_000000_V180529.cef")




