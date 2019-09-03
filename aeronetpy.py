# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 18:35:06 2018

@author: Menk
"""

#----------------------------------------------------------------------------

def convertDate(dates, inputformat = 'dd:mm:yyyy', outputformat = 'yyyy-mm-dd',JDNon = False):
    '''
    Function to convert date varible into user specified formate, also convert the 
    Gregorian calendar data to Julian Day Number (JDN)
    Input:
        dates, list like varible, each item in the list is a string
        imputformat, formate of the date varible
        outputformat, user specified date formate
        jdn, switch of Julian Day Number
    Output:
        date_output, output of user specified date formate
        jnd_output, Julian Day Number
    '''
    index_day = inputformat.index('dd')
    index_month = inputformat.index('mm')
    index_year = inputformat.index('yyyy')
    
    
    date_output = []
    jnd_output = []
    for date in dates:
        
        day = date[index_day:index_day+2]
        month = date[index_month:index_month+2]
        year = date[index_year:index_year+4]
        
        output_temp = outputformat
        output_temp = output_temp.replace('dd',day)
        output_temp = output_temp.replace('mm',month)    
        output_temp = output_temp.replace('yyyy',year)
        date_output.append(output_temp)
        if JDNon == True:
            jdn_temp = convertJDN(year+month+day)
            jnd_output.append(jdn_temp)
    if JDNon == True:
        return date_output, jnd_output
    else:
        return date_output

#----------------------------------------------------------------------------
def convertTime(dates, times, inputformat = 'hh:mm:ss'):
    '''
    Function to convert time varible string format into a float format
    Input:
        times, list like varible, each item in the list is a string
        imputformat, formate of the Time varible
    Output:
        float format time varible
    '''
    import datetime
    import numpy as np
    time_out = []
    timeLocal = []
    for date, time in zip(dates,times):      
    
        utcTime = datetime.datetime.strptime(str(date)+' '+str(time), "%d:%m:%Y %H:%M:%S")            
        time_float = (utcTime.hour*60.0*60.0 + utcTime.minute*60.0 + utcTime.second)/24.0/60.0/60.0
        time_out.append(time_float)

    time_out = np.array(time_out)
    return time_out

#----------------------------------------------------------------------------
def convertJDN(GregDay, outtype = 'int'):
    '''
    Function to convert the Gregorian calendar data to Julian Day Number (JDN)
    input format: yyyymmdd
    output: Julian Day Number, int or string, defalt type int
    '''
    import numpy as np
    year = np.int(GregDay[0:4])
    month = np.int(GregDay[4:6])
    day = np.int(GregDay[6:8])
    Julian_a = (14-month)//12
    Julian_y = year + 4800 - Julian_a
    Julian_m = month + 12 * Julian_a - 3
    
    jdn = day + (153*Julian_m+2)//5 + 365*Julian_y + Julian_y//4 - Julian_y//100 + Julian_y//400 -32045
    if outtype == 'str':
        jdn = str(jdn)        
    return jdn

#----------------------------------------------------------------------------
def convertGergDay(jnd,outputformat = 'yyyy-mm-dd'):
    '''
    http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    '''
    l = jnd + 68569
    n = 4*l//146097
    
    l= l-(146097*n+3)//4
    year= 4000*(l+1)//1461001
    l= l-1461*year//4+31
    month= 80*l//2447
    day = l-2447*month//80
    l= month//11
    month= month + 2 - 12*l
    year= 100*(n-49) + year + l   
    
    if day < 10:
        day = '0' + str(day)
    if month <10:
        month = '0' + str(month)
    GregDay = outputformat
    GregDay = GregDay.replace('dd',str(day))
    GregDay = GregDay.replace('mm',str(month))   
    GregDay = GregDay.replace('yyyy',str(year))   
    
    return GregDay

#----------------------------------------------------------------------------
def readAeronet(filename, vars, span = [], JDNon = True, outputformat = 'yyyy-mm-dd', keyword = 'Date(dd:mm:yyyy)' ):
	'''
	Fucntion to read the Aeronet files
	Input:
		filename, string like varible, specifiy file to be read
		var, list like varible, 
		span
		jdn
	Output:
		output, a dictionary like varible

	'''   
	import numpy as np
	import pandas as pd

	#use the Date(dd:mm:yyyy) item to detect the header line
	print('  - Readig: ', filename)
	headernames, count_header =  exploreAeronet(filename, keyword)
	data = pd.read_csv(filename, header = count_header, na_values =[-999])


	if len(span) == 1:     
		index = data.index[data[keyword] == span[0]].tolist()
		if len(index) != 0:   
			data = data.iloc[index,:]
		else:
			print('Required date has no data.')
		
	if len(span) == 2:
		index_s = data.index[data[keyword] == span[0]].tolist()
		index_e = data.index[data[keyword] == span[1]].tolist()       
		if (len(index_s) != 0) & (len(index_e) != 0):
			data = data.iloc[index_s[0]:index_e[-1],:]
		else:
			print('Required date has no data.')

	items = [keyword,'Time(hh:mm:ss)']
	items.extend(vars) 

	warningVar = ''
	for var in vars:
		if var not in headernames:
			warningVar = warningVar + ' '  + var
	if len(warningVar) > 0 :
		print('  Warning: ' + warningVar + ' is(are) not in the list...')
		
		return headernames


	exploreAeronet(filename, keyword = 'Date(dd:mm:yyyy)')
 
	output = {}   

	for item in items:   
		if item in headernames:            
			if item == keyword:
				date_out, jdn = convertDate(data[item].values, outputformat = outputformat, JDNon = JDNon)
				output ['Date' ] = np.array(date_out)
				output ['JDN' ] = np.array(jdn)
			elif item =='Time(hh:mm:ss)':
				output['TimeFloat'] = convertTime(data[keyword].values, data[item].values)

				output['Time'] = np.array((data[item].values))
			else:
				output[item] = np.array(data[item].values)
		else:
			print('Item', item, 'is not in the dataset.')    


	output['JDNLocal'] = output['JDN'] + np.array(output['TimeFloat'])

	jdnLocal = output['JDN'] + np.array(output['TimeFloat'])

	jdn_max = max(jdn)
	jdn_min = min(jdn)
	jdns = np.arange(jdn_min, jdn_max +1, 1,dtype = int)
	year = convertGergDay(jdn_min)[0:4]
	jdn_ini = convertJDN(year + '0101')

	output['JDNLocal'] = output['JDNLocal'] - jdn_ini + 1

	daily_index = []
	for junliandaynumber in jdns:
		index = np.where((jdnLocal >= junliandaynumber) & (jdnLocal <= junliandaynumber+1))[0]
		daily_index.append((junliandaynumber-jdn_ini+1,junliandaynumber,index))
	
	output['DailyIndex'] = daily_index  

	return output

#----------------------------------------------------------------------------
def exploreAeronet(filename, keyword = 'Date(dd:mm:yyyy)'): 
	import numpy as np 
	import pandas as pd
    #use the Date(dd:mm:yyyy) item to detect the header line
	with  open(filename,'r') as f:
		count_header = 0
		data = f.readline()
		while keyword not in data:
			data = f.readline()
			count_header = count_header + 1
	#read file
	data = pd.read_csv(filename, header = count_header, na_values =[-999])
	 
	validItem = []
	for item in data:
		temp_value = data[item].values
		count = np.where(temp_value == temp_value)
		if np.size(count) > 0:
			validItem.append(item)
	
	headernames = list(data.columns.values)  
	return validItem, count_header

#----------------------------------------------------------------------------
def selectData(sourcefile, parameters, refrence_parameter, aodThredhold, outputformat):
    
    import numpy as np
    import warnings
    warnings.filterwarnings("ignore",category =RuntimeWarning)
    
    jdn_s = sourcefile['DailyIndex'][0][1]
    jdn_e = sourcefile['DailyIndex'][-1][1]
    jdns = np.arange(0,jdn_e-jdn_s+1, dtype = int)
    

    keys = []
    for parameter in parameters:
        paired_parameter = 'DailyMean_' + parameter
        keys.append((paired_parameter, parameter))    
    
    print ('------------------')   
    Message = 'Using ' + refrence_parameter + ' as refrence.'
    print(Message)
    
    refrence_parameter = 'DailyMean_' + refrence_parameter
    
    #Add the new key in the dictionary
    for key in keys:
        sourcefile[key[0]]= []
    
    for jdn in jdns:
        index = sourcefile['DailyIndex'][jdn][2]
        for key in keys:
            sourcefile[key[0]].append(np.nanmean(sourcefile[key[1]][index]))

    refrence = np.array(sourcefile[refrence_parameter])
    
    refrence_indices = np.where(refrence <= aodThredhold)[0]
    diff = np.diff(refrence_indices)
    index_diff = np.where(diff == 1)[0]
    cleandayindex = refrence_indices[index_diff]
    
    var = parameters.copy()
    
    #add the default keyword in to the output dictionary
    var.append('Time')
    var.append('TimeFloat')
    output = {'Date':[],'JDN':[]}
    for item in var:
        output[item]= []

    for index in cleandayindex:
        jdn_1 = sourcefile['DailyIndex'][index][1]
        jdn_2 = jdn_1 + 1
        date_1 = convertGergDay(jdn_1,outputformat = outputformat)
        date_2 = convertGergDay(jdn_2,outputformat = outputformat)
        output['Date'].append((date_1,date_2))
        output['JDN'].append((jdn_1,jdn_2))
        index_1 = sourcefile['DailyIndex'][index][2]
        index_2 = sourcefile['DailyIndex'][index+1][2]
    
        for item in var:        
            aod_temp = (np.array(sourcefile[item])[index_1], np.array(sourcefile[item])[index_2])
            output[item].append(aod_temp)    
    
    return output, sourcefile 

#----------------------------------------------------------------------------
def writeCSV(dicData,filename):
    import pandas as pd      
    df = pd.DataFrame.from_dict(dicData,'columns')
    df.to_csv(filename, index=False)
    return df
    
#----------------------------------------------------------------------------
def get_tau(alpha, tau0, lamda0, lamda):
	
	tau = tau0 * (lamda / lamda0)**(-alpha) 
	
	return tau
