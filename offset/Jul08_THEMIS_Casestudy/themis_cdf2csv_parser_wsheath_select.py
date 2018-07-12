# -*- coding: utf-8 -*-
"""
Created on Fri Jun  8 13:39:40 2018

@author: Choon
"""


import os
import io
import matplotlib
import matplotlib.pyplot as plt
#import time
from datetime import datetime

#SETTING='fgmdata'
SETTING='spinaxis'

WRITING_SHEATH_INTERVALS=False

#cef_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
#cef_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"
#cef_file_name = "C3_CP_AUX_SPIN_AXIS__20060301_000000_20060302_000000_V130205"
#cef_file_name = "C3_CP_AUX_POSGSE_1M__20060301_000000_20060302_000000_V091201"
#cef_file_name = "C3_CP_FGM_5VPS__20061103_000000_20061104_000000_V140305"
#cef_file_name = "C3_CP_AUX_SPIN_AXIS__20060307_000000_20060312_000000_V130205"
#cef_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060312_000000_V140305"
#cef_file_name = "C3_CP_FGM_5VPS__20060312_000000_20060322_000000_V140305"
#cef_file_name = "C3_CP_FGM_5VPS__20060322_000000_20060401_010000_V140305"
cef_file_name = "C3_CP_AUX_SPIN_AXIS__20060322_000000_20060401_010000_V130205"

yyyy_s = '2006'
mm_s = '03'
dd_s = '30'
hh24_s = '15'

yyyy_e = '2006'
mm_e = '03'
dd_e = '30'
hh24_e = '20'


sheath_start_time = matplotlib.dates.date2num(datetime.strptime(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
sheath_end_time = matplotlib.dates.date2num(datetime.strptime(yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

sheathfile_name = "C3_CP_sheathintervals__200603"

if SETTING=='fgmdata':
    outputfile_name = "C3_CP_FGM_5VPS__200603_sheathselected"
elif SETTING=='spinaxis':
    outputfile_name = "C3_CP_AUX_SPIN_AXIS__200603_sheathselected"    

inputfile = open(os.getcwd()+"\\"+cef_file_name+".cef")
text = inputfile.readlines()
inputfile.close()


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
        if ((text[i]).split())[0] == 'DATA_UNTIL':
            break
    i+=1
    
data_start_pos = i+1

#########################################
#csv_file = io.open(os.getcwd() +"\\"+ file_name + ".csv",'a')
csv_file = io.open(os.getcwd() +"\\"+ outputfile_name + ".csv",'a')
if WRITING_SHEATH_INTERVALS:
    sheath_csv_file = io.open(os.getcwd() +"\\"+ sheathfile_name + ".csv",'a')

    #data_arr=[]

    sheath_csv_file.write(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z' +','+ yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z' +'\n')


i=data_start_pos
#i=3400000
while ((text[i]).split())[0] != '!RECORDS=':
    if SETTING == 'fgmdata':
        data_row_str = ((text[i]).split())[0]
        row_date = (data_row_str.split(','))[0]
        row_date_num = matplotlib.dates.date2num(datetime.strptime(row_date,'%Y-%m-%dT%H:%M:%S.%fZ'))
        if row_date_num > sheath_start_time and row_date_num < sheath_end_time:
            csv_file.write(data_row_str+'\n')
        if row_date_num > sheath_end_time:
            csv_file.write('sheath_interval_end'+'\n')
            print(i)
            break
    elif SETTING == 'spinaxis':
        data_row_str = (text[i])
        data_row_str.replace(" ","")
        #print(data_row_str)
        csv_file.write(data_row_str) 
        

    #print(text[i])
    
    i+=1








    

#ceflib.read(os.getcwd()+"\\C1_CP_FGM_SPIN__20180110_000000_20180117_000000_V180529.cef")




