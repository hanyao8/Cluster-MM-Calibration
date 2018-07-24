README


PREPROCESSING SETTINGS:

cef2csv:
setting: (fgmdata, spinaxis, ion)
cef_file_name: (The name of the cef file which is to be translated to csv.
					The file is to be placed in the same directory as 
					Preprocessing_Settings.txt and cef2csv_parser.py)

cef2csv_parser_wsheath_select:
setting: (fgmdata, spinaxis)
writing_sheath_intervals: (True,False) Set to 'True' to write sheath-pass
							and end times to the file to be named below.
cef_file_name: (The name of the cef file which is to be translated to csv.
					The file is to be placed in the same directory as 
					Preprocessing_Settings.txt and cef2csv_parser.py)
yyyy_s = Single sheath-pass starting time- year (yyyy)
mm_s = Single sheath-pass starting time- month (mm)
dd_s = Single sheath-pass starting time- day (dd)
hh24_s = Single sheath-pass starting time- hour (hh 24 e.g. 23 for 11 in evening)
yyyy_e = Single sheath-pass ending time- year (yyyy)
mm_e = Single sheath-pass ending time- month (mm)
dd_e = Single sheath-pass ending time- day (dd)
hh24_e = Single sheath-pass ending time- hour (hh 24)
sheathfile_name = File to place sheath-pass starting and ending times
outputfile_name = File to place sheath selected FGM data

gse2dsl_transformation
themis_validation: (True,False) Set to 'False' for typical analysis of Cluster data
gsedata_file_name: Name of GSE-based FGM data file to be transformed
spinaxis_file_name: Name of spin-axis orientation time series file to be used for transformation


ANALYSIS SETTINGS:
ARTIFICIAL_OFFSET: (True,False)
PLOT: (True,False) 'True' to generate graphs
artificial_Bz_offset: (Float) This artificial offset is applied if ARTIFICIAL_OFFSET==True
csv_file_name: DSL based FGM data file
sheath_file_name: Name of file containing the sheath-pass starting and ending times
yyyy_s = Overall interval starting time- year (yyyy)
mm_s = Overall interval starting time- month (mm)
dd_s = Overall interval starting time- day (dd)
hh24_s = Overall interval starting time- hour (hh 24 e.g. 23 for 11 in evening)
yyyy_e = Overall interval ending time- year (yyyy)
mm_e = Overall interval ending time- month (mm)
dd_e = Overall interval ending time- day (dd)

