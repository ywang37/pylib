'''
# package of download data from NASA LAADS...

example: 
from pylib.dataDownloadpy import *

rootDir = '/Users/mzhou16/project/processes_annual_fire/pythonCode/'
products = ['VNP14_N', 'VNP03DNB_N']
date = ['2019-03-01', '2019-03-1']
cord = [50, 40, -103.8, -102.2]
download_data(products, date, cord, rootDir, crosscheck = True)

'''
from mylib.merra2.download import *
#-----------------------------------------------------------------------
def	create_dir(rootDir, dirDic):
	import os
	for key in dirDic.keys():
		subDir = dirDic[key][0]
		dir = rootDir + key + '/'
		if os.path.isdir(dir) == False:
			os.mkdir(dir)
			print('  - Making run directory: ' + dir)	
	return

#-----------------------------------------------------------------------
def get_data(product, date, cord, rootDir, **kwargs):

	import xml.dom.minidom
	import pandas as pd
	import requests
	import os
	import time
	import json
	
	fileNum = 0
	obtainedFileNum = 0
	
	print(' - get_data:', product)
	
	appKey	= kwargs.get('appKey', 'Bearer B6148868-1922-11E9-A62B-CCB570C49BBF')
	
	outputDir = rootDir + product[0] + '/'
	
	# create the time string 
	startTime = date[0]
	endTime   = date[1]
	period = pd.date_range(start=startTime, end=endTime, freq='D')
	period_string = period.strftime('%Y-%m-%d')
	
	filenames = {}
	
	productname = product[0]
	
	for item_date in period_string:
	
		urls = get_data_url(product, item_date, cord)
		
		for url in urls:
			filename = url.split('/')[-1]
			filenames[filename] = url
			#-------------------
			# core download...
			download_singleFile(url, outputDir, **kwargs)
			
	outfilename = rootDir + productname + '.json'
	with open(outfilename, 'a') as outfile:
		json.dump(filenames, outfile)
	
	fileNum = len(filenames)
	
	obtainedFilenames = os.listdir(outputDir)
	obtainedFileNum = len(obtainedFilenames)
	
	# we do empty file and missing file check...
	while (fileNum != obtainedFileNum):
		print( ' - numer of file != obtained numer of file')
		obtainedFilenames = os.listdir(outputDir)
		obtainedFileNum = len(obtainedFilenames)
		
		count = 0
		
		for filename in filenames.keys():
		
			if filename not in obtainedFilenames:
				print(filename)
				print(' - Missed file: ', filename)
				download_singleFile(filenames[filename], outputDir, **kwargs)
			else:

				filesize = os.path.getsize(outputDir + filename)
				countEmpt = 0
							
				while (filesize < 100):
					print(' - Empty file: ', filename)
					download_singleFile(filenames[filename], outputDir, **kwargs)
					filesize = os.path.getsize(outputDir + filename)
					countEmpt = countEmpt + 1
					if (countEmpt > 5):
						print(' - get_data: Warning! Download process is aborted by empty data')
						break
		
		count = count + 1
		if (count > 0):
			print(' - get_data: Warning! Download process is aborted by empty data')
			print(filename)
			break
			
	print(' - END OF DOWNLOADING!!!')			
	print(' - Total number of files that need to be downloaded: ', fileNum)
	print(' - Total numer of downloaded file: ', obtainedFileNum)
	
	return fileNum, obtainedFileNum

#-----------------------------------------------------------------------
def get_satDic(products):

	'''
	Python code to return a product dictionary
	
	Parameters:
	products: list, name of the product, refer to the key of the dirDic below
	
	Return:
	outputDic: product dictionary, contains the information of the download data	
	https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/
	
	'''

	dirDic = {}
	
	#------------------
	# VIIRS DNB
	#------------------
	dirDic['VNP02DNB_N'] = ['VNP02DNB_N', 'VNP02DNB', '5110', 'N']
	dirDic['VNP03DNB_N'] = ['VNP03DNB_N', 'VNP03DNB', '5110', 'N']
	
	dirDic['VNP02DNB_D'] = ['VNP02DNB_D', 'VNP02DNB', '5110', 'D']
	dirDic['VNP03DNB_D'] = ['VNP03DNB_D', 'VNP03DNB', '5110', 'D']
	
	dirDic['VNP02DNB_B'] = ['VNP02DNB_N', 'VNP02DNB', '5110', 'NB']
	dirDic['VNP03DNB_B'] = ['VNP03DNB_N', 'VNP03DNB', '5110', 'NB']

	#------------------
	# VIIRS M band
	#------------------
	dirDic['VNP02MOD_N'] = ['VNP02MOD_N', 'VNP02MOD', '5110', 'N']
	dirDic['VNP03MOD_N'] = ['VNP03MOD_N', 'VNP03MOD', '5110', 'N']

	dirDic['VNP02MOD_D'] = ['VNP02MOD_D', 'VNP02MOD','5110', 'D']
	dirDic['VNP03MOD_D'] = ['VNP03MOD_D', 'VNP03MOD','5110', 'D']
	
	dirDic['VNP02MOD_B'] = ['VNP02MOD_N', 'VNP02MOD','5110', 'NB']
	dirDic['VNP03MOD_B'] = ['VNP03MOD_N', 'VNP03MOD','5110', 'NB']

	#------------------
	# VIIRS I band
	#------------------

	dirDic['VNP02IMG_N'] 	= ['VNP02IMG_N', 'VNP02IMG', '5110', 'N']
	dirDic['VNP03IMG_N'] 	= ['VNP03IMG_N', 'VNP03IMG', '5110', 'N']

	dirDic['VNP02IMG_B'] 	= ['VNP02IMG_N', 'VNP02IMG', '5110', 'NB']
	dirDic['VNP03IMG_B'] 	= ['VNP03IMG_N', 'VNP03IMG', '5110', 'NB']
	
	#------------------------
	# VIIRS activate fire
	#------------------------

	dirDic['VNP14_N']	 = ['VNP14_N', 'VNP14','5000', 'N']

	dirDic['VNP14_D'] 	 = ['VNP14_D', 'VNP14','5000', 'D']
	
	dirDic['VNP14_B'] 	 = ['VNP14_N', 'VNP14','5000', 'NB']

	dirDic['VNP14IMG_N'] 	= ['VNP14IMG_N', 'VNP14IMG', '5000', 'N']
	
	dirDic['VNP14IMG_D'] 	= ['VNP14IMG_D', 'VNP14IMG', '5000', 'D']

	dirDic['VNP14IMG_B'] 	= ['VNP14IMG_N', 'VNP14IMG', '5000', 'B']
	
	#------------------------------------
	# VIIRS black marble nighttime light
	#------------------------------------
	
	dirDic['VNP46A1']    = ['VNP46A1', 'VNP46A1', '5000', 'N']
	
	#---------------------------
	# MODIS surface reflectance
	#---------------------------
	
	dirDic['MCD43C3']           = ['MCD43C3', 'MCD43C3', '6', 'D']
	
	#-----------------
	# MODIS MAIAC AOD
	#-----------------
	dirDic['MCD19A2'] 		= ['MCD19A2',  'MCD19A2',   '6', 'D']
	dirDic['MOD02HKM'] 		= ['MOD02HKM', 'MOD02HKM', '61', 'D']
	dirDic['MOD021KM'] 		= ['MOD021KM', 'MOD021KM', '61', 'D']
	dirDic['MOD03'] 		   = ['MOD03',    'MOD03',    '61', 'D']
	dirDic['MOD04_3K'] 		= ['MOD04_3K', 'MOD04_3K', '61', 'D']

	keys = dirDic.keys()
	outputDic = {}
	for product in products:
		if product in keys:
			outputDic[product] = dirDic[product]
	
	return outputDic

#-----------------------------------------------------------------------
def get_data_url(product, date, cord):

	import xml.dom.minidom
	import pandas as pd
	import requests

	north = str(cord[0])
	south = str(cord[1])
	west  = str(cord[2])
	east  = str(cord[3])

	productName = product[1]
	collection  = product[2]
	DayNigStatus= product[3]
	urls=[]
	url='https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/searchForFiles' + \
		'?product=' + productName + \
		'&collection=' + collection + \
		'&start=' + date + \
		'&stop='  + date + \
		'&north=' + north     + \
		'&south=' + south     + \
		'&west='  + west      + \
		'&east='  + east      + \
		'&coordsOrTiles=' + 'coords'+ \
		'&dayNightBoth='      + DayNigStatus
	print(' - Sending request to',url)
	response = requests.get(url) 
	doc = xml.dom.minidom.parseString(response.content)
	response_fields=[]
	for node in doc.getElementsByTagName("return"):
		response_fields.append(node.firstChild.nodeValue)

		
	fileIds = ",".join(response_fields)
	fileIds = fileIds#.encode('utf-8')
	url='https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/getFileUrls?fileIds='+ \
		 fileIds
	response = requests.get(url)
	doc = xml.dom.minidom.parseString(response.content)
	
	for node in doc.getElementsByTagName("return"):

		urls.append(node.firstChild.nodeValue)
							
	return urls
	
#-----------------------------------------------------------------------
def download_singleFile(item, outputDir, **kwargs):

	import os
	import time
	import numpy as np
	from random import randint
	
	appKey	= kwargs.get('appKey', 'Bearer B6148868-1922-11E9-A62B-CCB570C49BBF')
	
	# -q we mute wget command to print any message to the screen
	command = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 ' + item + \
			  ' --header "Authorization: ' + appKey + '" -O ' + \
			  outputDir + os.path.basename(item)

	print (' - Downloading', item)
	os.system(command)
	numSecond = np.float(randint(5,10))/10
	time.sleep(numSecond)
	
	count = 0	
	filesize = os.path.getsize(outputDir + os.path.basename(item))
	while filesize < 100:
		os.system(command)
		
		filesize = os.path.getsize(outputDir + os.path.basename(item))
		
		count = count + 1
		if count>5:
			break
	return

#-----------------------------------------------------------------------
def download_data(products, date, cord, rootDir, **kwargs):
	
	appKey		= kwargs.get('appKey', 'Bearer B6148868-1922-11E9-A62B-CCB570C49BBF')
	crosscheck	= kwargs.get('crosscheck', False)
	
	dirDic = get_satDic(products)
	
	create_dir(rootDir, dirDic)
	
	fileNums 		= [] 
	obtainedFileNums = []
	productNames 	= []
	
	for key in dirDic.keys():
		fileNum, obtainedFileNum = get_data(dirDic[key], date, cord, rootDir, **kwargs)
		
		fileNums.append(fileNum)
		obtainedFileNums.append(obtainedFileNum)
		productNames.append(key)
		
	print()
	print(' - download_data summary')
	for fileNum, obtainedFileNum, productName in zip(fileNums, \
	                                                 obtainedFileNums, \
	                                                 productNames):
		print(' - product name:', productName, \
		      ', file number: %i, actual download number: %i'\
		      %(fileNum, obtainedFileNum))
	
	if crosscheck:
		missfilenames  = download_missing(products, rootDir, **kwargs)
		emptyfilenames = download_empty(products, rootDir, **kwargs)	
	
#-----------------------------------------------------------------------
def download_missing(products, rootDir, **kwargs):
	
	import os
	import mylib.convert_time_coordinate as ctc
	
	appKey		= kwargs.get('appKey', 'Bearer B6148868-1922-11E9-A62B-CCB570C49BBF')
	
	numDir = len(products)
	allFiles = {}
	
	dirDic = get_satDic(products)

	for dirName in products:
		filenames = os.listdir(rootDir + dirName + '/')
		print(' - download_missing:', products, len(filenames))
		for filename in filenames:
			if '.nc' not in filename:
				continue
			tempSTR = filename.split('.')
			product = dirName #tempSTR[0]
			jdn		= tempSTR[1]
			orbit	= tempSTR[2]
			# we use the jdn and orbit and key...
			key = (jdn, orbit)	
			if key not in allFiles.keys():
				allFiles[key] = []
				allFiles[key].append(product)
			else:
				allFiles[key].append(product)
	
	
	missfilenames = []
	cord = [90, -90, -180, 180]
	for key in allFiles.keys():
		# find the each day and each orbit and count the number of the files
		# if the number is less than the number of directories, find out
		# which file is missing...
		if len(allFiles[key]) < numDir:
			
			for product in products:
				
				if product not in allFiles[key]:
					# get the day and orbit of the file
					
					jdn = key[0]
					orbit = key[1]
					date = ctc.get_displatDate(jdn, format = 'yyyy-mm-dd')
					
					urls = get_data_url(dirDic[product], date, cord)
					
					# set the identity of the file...
					strID = jdn + '.' + orbit
				
					# download the file
					for item in urls:
						
						if strID in item:
							saveDir = rootDir + product + '/'
							print( ' - The missing file is', item)
							download_singleFile(item, saveDir, **kwargs)
							missfilenames.append(item)
	print()
	if len(missfilenames) == 0:
		print(' - download_missing: No missing files')
	else:
		print(' - download_missing: number of missing file downloaded', len(missfilenames))
	
	return missfilenames

#-----------------------------------------------------------------------
def download_empty(products, rootDir, **kwargs):

	import os
	appKey		= kwargs.get('appKey', 'Bearer B6148868-1922-11E9-A62B-CCB570C49BBF')
	numDir = len(products)
	dirDic = get_satDic(products)
	urlhead = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/'

	emptyFile = {}
	
	for dirName in products:
		emptyFile[dirName] = []

		filenames = os.listdir(rootDir + dirName +'/')
		for filename in filenames:
			size = os.path.getsize(rootDir + dirName +'/' + filename)
			if size < 400:
				emptyFile[dirName].append(filename)
	
	emptyfilenames = []
	for key in emptyFile.keys():
		for filename in emptyFile[key]:
			print( ' - The empty file is', filename)
			tempSTR = filename.split('.')
			year    = tempSTR[1][1:5]
			jdn		= tempSTR[1][5:]
			
			productName = dirDic[key][1]
			collection  = dirDic[key][2]
			
			webDir  = urlhead + collection + '/' + productName + '/' + \
			          year + '/' + jdn + '/'

			item    = webDir + filename	
			outputDir = rootDir + key + '/'
			download_singleFile(item, outputDir, **kwargs)
			emptyfilenames.append(item)

	print()
	if len(emptyfilenames) == 0:
		print(' - download_empty: No empty files')
	else:
		print(' - download_empty: number of empty file downloaded', len(emptyfilenames))
		
	return emptyfilenames

#-----------------------------------------------------------------------
def get_data_properties(product, date, cord, rootDir):

	import xml.dom.minidom
	import pandas as pd
	import requests
	import xml.etree.ElementTree as ET
	import os
	import time
	
	
	bad_file = []
	cont_lim = 2
	outputDir = rootDir + product[0] + '/'
	
	# create an directory for the data...
	if os.path.isdir(outputDir) == False:
		os.mkdir(outputDir)
		print('  - Making run directory: ' + outputDir)	
	
	#
	startTime = date[0]
	endTime   = date[1]
	period = pd.date_range(start=startTime, end=endTime, freq='D')
	period_string = period.strftime('%Y-%m-%d')

	north = str(cord[0])
	south = str(cord[1])
	west  = str(cord[2])
	east  = str(cord[3])

	productName = product[1]
	collection  = product[2]
	DayNigStatus= product[3]


	for item_date in period_string:
		url='https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/searchForFiles' + \
			'?product=' + productName + \
			'&collection=' + collection + \
			'&start=' + item_date + \
			'&stop='  + item_date + \
			'&north=' + north     + \
			'&south=' + south     + \
			'&west='  + west      + \
			'&east='  + east      + \
			'&coordsOrTiles=' + 'coords'+ \
			'&dayNightBoth='      + DayNigStatus
		print('  - Sending request to', url)
		response = requests.get(url)  
		doc = xml.dom.minidom.parseString(response.content)
		response_fields=[]
		for node in doc.getElementsByTagName("return"):
			response_fields.append(node.firstChild.nodeValue)
		fileIds = ",".join(response_fields)
		fileIds = fileIds#.encode('utf-8')
		
		url='https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/getFileUrls?fileIds='+ \
		     fileIds
		response = requests.get(url)
		doc = xml.dom.minidom.parseString(response.content)
		urls=[]
		for node in doc.getElementsByTagName("return"):
			urls.append(node.firstChild.nodeValue)
		
		urlProp='https://modwebsrv.modaps.eosdis.nasa.gov/axis2/services/MODAPSservices/getFileProperties?fileIds='+fileIds
		responseProp=requests.get(urlProp)
		
		tree = ET.ElementTree(ET.fromstring(responseProp.text))
		root = tree.getroot()
		
		filesInfo = []
		for child in root:
			fileInfo = {'checksum': None , 'fileId': None, 'fileName': None, 'fileSizeBytes': None, 
				 		'fileType': None, 'ingestTime': None, 'online': None, 'startTime': None, 'urlFile': None}
			for child1 in child:
				field = child1.tag.replace('{http://modapsws.gsfc.nasa.gov/xsd}', '')
				fileInfo[field] = child1.text
				if field == 'fileName':
					url_aux = [ url for url in urls if child1.text in url]
					if len(url_aux) == 1:
						fileInfo['urlFile']= url_aux[0]	
			filesInfo.append(fileInfo)

		for aux_info in filesInfo:
			if aux_info['online'] == 'true' and int(aux_info['fileSizeBytes']) > 1800 :
				bnameF = os.path.basename(aux_info['urlFile'])
				if os.path.isfile(outputDir + bnameF) and (os.path.getsize(outputDir + bnameF) == int(aux_info['fileSizeBytes'])):
					# the file is already in the output directory
					print(' - File already exited...SKIP!!')
					continue
				
				command = 'wget -e robots=off -m -np -R .html,.tmp -nH --cut-dirs=3 ' + aux_info['urlFile'] + \
			          ' --header "Authorization: Bearer B6148868-1922-11E9-A62B-CCB570C49BBF" -O ' + \
			          outputDir + bnameF #+ " > /dev/null 2>&1"
			    			
				os.system(command)
				
				counter = 0
				while ( os.path.getsize(outputDir + bnameF) != int(aux_info['fileSizeBytes'])):
					print( ' - Find broken file, redownloading...')
					os.system(command)
					counter = counter + 1
					if counter>cont_lim:
						bad_file.append(aux_info['urlFile'])
						break
				time.sleep(0.4)

	if len(bad_file)>0:
		with open(str(time.time()).split('.') + 'bad_download.txt', 'w') as f:
			for item in bad_file:
				f.write("%s\n" % item)

					
	return








