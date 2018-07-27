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
import numpy as np
from pathlib import Path

preprocessing_settings = open(os.getcwd()+"\\"+"Preprocessing_Settings"+".txt")
pp_text = preprocessing_settings.readlines()

#SETTING='fgmdata'
#SETTING='spinaxis'

if (pp_text[6].split())[0] != 'cef2csv_parser_wsheath_select':
    raise Exception("Error: Check Preprocessing_Settings.txt formatting")

SETTING = ( (pp_text[7]).split() )[1]

#WRITING_SHEATH_INTERVALS=True
WRITING_SHEATH_INTERVALS= bool((pp_text[8].split())[1])

#cef_file_name = "C1_CP_FGM_5VPS__20060301_103000_20060301_113000_V140304"
#cef_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060302_000000_V140305"
#cef_file_name = "C3_CP_AUX_SPIN_AXIS__20060301_000000_20060302_000000_V130205"
#cef_file_name = "C3_CP_AUX_POSGSE_1M__20060301_000000_20060302_000000_V091201"
#cef_file_name = "C3_CP_FGM_5VPS__20061103_000000_20061104_000000_V140305"
#cef_file_name = "C3_CP_AUX_SPIN_AXIS__20060307_000000_20060312_000000_V130205"
#cef_file_name = "C3_CP_FGM_5VPS__20060301_000000_20060312_000000_V140305"
#cef_file_name = "C3_CP_FGM_5VPS__20060312_000000_20060322_000000_V140305"
#cef_file_name = "C3_CP_FGM_5VPS__20060322_000000_20060401_010000_V140305"
#cef_file_name = "C3_CP_FGM_5VPS__20051021_000000_20051031_230000_V140305"
cef_file_name = ( (pp_text[9]).split() )[1]

#yyyy_s = '2005'
#mm_s = '10'
#dd_s = '21'
#hh24_s = '00'

#yyyy_e = '2005'
#mm_e = '10'
#dd_e = '31'
#hh24_e = '23'


yyyy_s = ( (pp_text[10]).split() )[1]
mm_s = ( (pp_text[11]).split() )[1]
dd_s = ( (pp_text[12]).split() )[1]
hh24_s = ( (pp_text[13]).split() )[1]

yyyy_e = ( (pp_text[14]).split() )[1]
mm_e = ( (pp_text[15]).split() )[1]
dd_e = ( (pp_text[16]).split() )[1]
hh24_e = ( (pp_text[17]).split() )[1]


sheath_start_time = matplotlib.dates.date2num(datetime.strptime(yyyy_s+'-'+mm_s+'-'+dd_s+'T'+hh24_s+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))
sheath_end_time = matplotlib.dates.date2num(datetime.strptime(yyyy_e+'-'+mm_e+'-'+dd_e+'T'+hh24_e+':00:00.000Z','%Y-%m-%dT%H:%M:%S.%fZ'))

if sheath_start_time > sheath_end_time:
    raise Exception("Error: Starting time of sheath interval is after the ending time")

sheathfile_name = ( (pp_text[18]).split() )[1]

if SETTING=='fgmdata':
    outputfile_name = ( (pp_text[19]).split() )[1]
elif SETTING=='spinaxis':
    outputfile_name = ( (pp_text[19]).split() )[1]

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


depfile_path = Path(os.getcwd() +"\\"+ "__temp_data_end_pos" + ".csv")

if depfile_path.is_file():
    print("1")
    depfile_r = open(os.getcwd() +"\\"+ "__temp_data_end_pos" + ".csv")
    depfile_text = depfile_r.readlines()  
    if depfile_text==[]:
        print("11")
        i=data_start_pos
    elif (depfile_text[-1].strip()).split(",")[0]==cef_file_name:
        print("12")
        last_end =int( (depfile_text[-1].strip()).split(",")[1] )
        i= int( last_end/(10**(int(np.log10(last_end)-1))) ) * (10**(int(np.log10(last_end)-1)))
    else:
        print("13")
        i=data_start_pos
else:
    print("2")
    i=data_start_pos

    
    





#i=3300000




print(i)
while ((text[i]).split())[0] != '!RECORDS=':
    if i%100000==0:
        print(i)
    if SETTING == 'fgmdata':
        data_row_str = ((text[i]).split())[0]
        row_date = (data_row_str.split(','))[0]
        row_date_num = matplotlib.dates.date2num(datetime.strptime(row_date,'%Y-%m-%dT%H:%M:%S.%fZ'))
        if row_date_num > sheath_start_time and row_date_num < sheath_end_time:
            csv_file.write(data_row_str+'\n')
        if row_date_num > sheath_end_time:
            #csv_file.write('sheath_interval_end'+'\n')
            print(i)
            data_end_pos = i
            break
    elif SETTING == 'spinaxis':
        data_row_str = (text[i])
        data_row_str.replace(" ","")
        #print(data_row_str)
        csv_file.write(data_row_str)
        

    #print(text[i])
    
    i+=1


depfile_w = io.open(os.getcwd() +"\\"+ "__temp_data_end_pos" + ".csv",'a')
depfile_w.write( cef_file_name +","+str( data_end_pos ) +"\n")



csv_file.close()
sheath_csv_file.close()
depfile_r.close()
depfile_w.close()



    

#ceflib.read(os.getcwd()+"\\C1_CP_FGM_SPIN__20180110_000000_20180117_000000_V180529.cef")




