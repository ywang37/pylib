#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-----------------------------------------------------------------------
def JulianDay(GregDay, **kwargs):
	'''
	Function to convert the Gregorian date to Julian Day Number(JDN)
	
	Reference: https://en.wikipedia.org/wiki/Julian_day

	Parameters
	----------
	GregDay	: string
			  format - yyyymmdd
			  
	Optional Parameters
	----------
	outtype : string, specify the output type, int, string, or nasa format
	type	: string
			  global - output is global Julian Day Number since November 23, −4713
			  local  - output is localized to that year
	Return
	----------    
	jdn		: int or str, defalt type int
			  Julian Day Number
	'''
	import numpy as np
	outtype = kwargs.get('outtype', 'int')
	type = kwargs.get('type', 'global')


	year = np.int(GregDay[0:4])
	month = np.int(GregDay[4:6])
	day = np.int(GregDay[6:8])
	Julian_a = (14-month)//12
	Julian_y = year + 4800 - Julian_a
	Julian_m = month + 12 * Julian_a - 3


	jdn = day + (153*Julian_m+2)//5 + 365*Julian_y + Julian_y//4 - Julian_y//100 + Julian_y//400 -32045


	if type == 'local':
		jdn_base = JulianDay( str(year) + '0101', type = 'global' )
		jdn = jdn - jdn_base + 1

	if outtype == 'str':
		jdn = np.str(jdn)

	if outtype == 'nasa':
		jdn_base = JulianDay( str(year) + '0101', type = 'global' )
		jdn = jdn - jdn_base + 1
		jdn = str(jdn)
		while len(jdn) < 3:
			jdn = '0' + jdn
		jdn = 'A' + GregDay[0:4] + jdn
	
	return jdn

#-----------------------------------------------------------------------
def GregorianDay(jnd, outputformat = 'yyyy-mm-dd', **kwargs):
    '''
    Function to convert the Julian Day Number(JDN)the Gregorian date
    http://aa.usno.navy.mil/faq/docs/JD_Formula.php

	Parameters
	----------
	jdn				: julian day of year
	outputformat	: number of day of the year

	Return
	----------    
	GregDay			: Gregorian date
	
    '''
    
    noYear = kwargs.get('noYear', False)
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
    if noYear:
    	GregDay = GregDay.replace('yyyy-','')
    else:
    	GregDay = GregDay.replace('yyyy',str(year))

    return GregDay

#-----------------------------------------------------------------------
def get_displatDate(date, format = 'mm-dd-yyyy', **kwargs):

	year = date[1:5]
	yearlyJDN = date[5:]	
	jdbBase = year + '0101'
	globalJDN = JulianDay(jdbBase,outtype = 'int') + int(yearlyJDN) - 1
	displatDate = GregorianDay(globalJDN, outputformat = format, **kwargs)
	
	return displatDate
	
#-----------------------------------------------------------------------
def get_date_series(gregDayBeg, gregDayEnd, **kwargs):

	import numpy as np
	outtype = kwargs.get('outtype', 'jdn')
	outputformat = kwargs.get('outputformat', 'yyyy-mm-dd')
	
	jdnBeg = JulianDay(gregDayBeg)
	jdnEnd = JulianDay(gregDayEnd) + 1
	jdns = np.arange(jdnBeg, jdnEnd)
	
	dateSeries = []
	if outtype == 'jdn':
		dateSeries = list( jdns )
		
	if outtype == 'greg':
		for jdn in jdns:
			GregDay = GregorianDay(jdn, outputformat = outputformat)
			dateSeries.append(GregDay)
			
	if outtype == 'nasa':
		for jdn in jdns:
			GregDay = GregorianDay(jdn, outputformat = 'yyyymmdd')
			dateSeries.append(JulianDay(GregDay, outtype = 'nasa'))		
	
	return dateSeries
#-----------------------------------------------------------------------
def get_month(dateBeg, dateEnd):

	import numpy as np
	monthDic = { 1  : 'JAN', 2  : 'FEB', 3  : 'MAR',
	             4  : 'APR', 5  : 'MAY', 6  : 'JUN',
	             7  : 'JUL', 8  : 'AUG', 9  : 'SEP',
	             10 : 'OCT', 11 : 'NOV', 12 : 'DEC'}
	            	
	monthBeg = dateBeg[4:6]
	monthEnd = dateEnd[4:6]
	
	monthBegInt = np.int(np.float(monthBeg))
	monthEndInt = np.int(np.float(monthEnd))
	
	monthStr = []
	
	for id in range(monthBegInt, monthEndInt + 1):
		monthStr.append(monthDic[id])
	
	return  monthStr

#-----------------------------------------------------------------------
def cal_grid(tile, numCeil, proj):
	'''
	Function to calculate the geographical coordinates for a given projection
	method
	
	Parameters
    ----------
		tile - str format, example: 'h07v05'
		numCeil - number of the ceils in one tile
		proj : projection method, Sinusoidal and PlateCarree available..
	Return
    ---------- 
		geographical coordinates of the tile
	'''
	if proj == 'Sinusoidal':
		latitude, longitude = cal_sinu_grid(tile, numCeil)
	if proj == 'PlateCarree':
		latitude, longitude = cal_PlateCarree_grid(tile, numCeil)
		
	return latitude, longitude
	
#-----------------------------------------------------------------------
def get_cord(cord, numCeil, proj):
	'''
	Function to get the geographical coordinates of the given region for 
	a given projection method

	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates 
			   (Top Latitude, Bottom Latitude, Left Longitude, Right Longitude)
		
		numCeil : number of pixel in one tile, default value 1200
		proj : projection method, Sinusoidal and PlateCarree available..
		
	Return
    ----------
		meshLat : array
				  longitude 
		meshLon : array
				  longitude 
		hidMin	: int
    		      vetical index of the point 
		hidMax	: int
    		  	 vetical index of the point 
		vidMin  : int
    		  	  vetical index of the point 
		vidMax  : int
    		  	  vetical index of the point
	'''
	if proj == 'Sinusoidal':
		meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = \
		get_cord_Sinusoidal(cord, numCeil, debug = 0)
	if proj == 'PlateCarree':
		meshLat, meshLon, hidMin, hidMax, vidMin, vidMax = \
		get_cord_PlateCarree(cord, numCeil, debug = 0)
		
	return meshLat, meshLon, hidMin, hidMax, vidMin, vidMax
	
#-----------------------------------------------------------------------
def get_point_tile(cord, pos):
	'''
	Function of calculating the tile of a specific point
	
	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates (latitude, longituede)
		
		pos : position of the point
		
	Return
    ----------
    	hid : int
    		  horizental index of the point
    	vid : int
    		  vetical index of the point
    	tile: str
    		  tile name
	'''
	if proj == 'Sinusoidal':
		latitude, longitude = get_point_tile_Sinusoidal(cord, pos)
	if proj == 'PlateCarree':
		latitude, longitude = get_point_tile_PlateCarree(cord, pos)
		
	return latitude, longitude


#-----------------------------------------------------------------------
# here after, function related to each projection method
'''
Function of Sinusoidal projection...
'''
def cal_sinu_grid(tile, numCeil):
	'''
	
	Function to calculate the geographical coordinates of the sinusoidal grid
	
	Parameters
    ----------
		tile - str format, example: 'h07v05'
		numCeil - number of the ceils in one tile
	
	Return
    ---------- 
		geographical coordinates of the tile
	
	Reference: 1. https://code.env.duke.edu/projects/mget/wiki/SinusoidalMODIS
			   2. https://onlinelibrary.wiley.com/doi/pdf/10.1111/0033-0124.00327
	'''
	import numpy as np

	halfCeilLen = 926.62543305/2.0
	halfHoriLenght = 20015109.354
	halfVertLenght = 10007554.677	
	numHoriTail = 37
	numVertTail = 19
	
	xx = np.linspace(-halfHoriLenght, halfHoriLenght, numHoriTail)
	yy = np.linspace(halfVertLenght, -halfVertLenght, numVertTail) 
		
	vid = np.int(np.float(tile[4:6]))
	hid = np.int(np.float(tile[1:3]))
	print('  - Calulating geographical coordinates of ', tile)
	print('    Vertical Tile:', vid, 'Horizontal Tile:', hid)

	x = np.linspace(xx[hid] + halfCeilLen, xx[hid+1] - halfCeilLen, numCeil)
	y = np.linspace(yy[vid] - halfCeilLen, yy[vid+1] + halfCeilLen, numCeil)
	
	xv, yv = np.meshgrid(x, y)
	
	latitude, longitude = sinu_to_geog((xv, yv))
	
	return latitude, longitude

#-----------------------------------------------------------------------
def sinu_to_geog(cord):
	'''
	Function of converting the sinusoidal projection to geographical projection
	
	Parameters
    ----------
		sinusoidal point - list or tuple liked, (x, y)
		
	Return
    ----------
		geographical coordinates - latitude, longituede

	'''
	import numpy as np
	x = cord[0]
	y = cord[1]
	pi = 180.0 / np.pi
	R = 6371007.181000
	
	phi = y/R
	lamda = x / np.cos(phi) / R
	
	latitude = phi * pi
	longituede = lamda * pi

	return latitude, longituede

#-----------------------------------------------------------------------
def geog_to_sinu(cord):
	'''
	Function of converting the geographical projection to sinusoidal projection
	
	Parameters
    ----------
		geographical coordinates - list or tuple liked, (latitude, longituede)
		
	Return
    ----------
		sinusoidal point - (x, y)
	
	
	'''
	
	import numpy as np

	lat = cord[0]
	lon = cord[1]
	
	pi = 180.0 / np.pi
	R = 6371007.181000
	
	phi = lat / pi
	lamda = lon / pi
	y = phi * R
	x = np.cos(phi) * lamda * R
	
	return x, y

#-----------------------------------------------------------------------
def get_point_tile_Sinusoidal(cord, pos):
	'''
	Function of calculating the tile of a specific point
	
	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates (latitude, longituede)
		
		pos : position of the point
		
	Return
    ----------
    	hid : int
    		  horizental index of the point
    	vid : int
    		  vetical index of the point
    	tile: str
    		  tile name
	'''
	import numpy as np
	

	halfHoriLenght = 20015109.354
	tileHoriLenght = halfHoriLenght/18
	
	halfVertLenght = 10007554.677
	tileVertLenght = halfVertLenght/9
	

	CeilLen = 926.625433056
	numCeil = 1200
	
	x, y = geog_to_sinu(cord)
	x = x + CeilLen/2.0
	y = y - CeilLen/2.0

	x_res = abs(np.round(x / tileHoriLenght) * tileHoriLenght - x)
	y_res = abs(np.round(y / tileVertLenght) * tileVertLenght - y)
	
	hid =  np.round(x // tileHoriLenght) + 18
	vid = 8 -  np.round(y // tileVertLenght)

	if y_res < CeilLen/2.0:
# 		print('yyyy')
		if pos == 'UpperLeft':
			vid = int(vid) + 1
		if pos == 'LowerRight':
			vid = int(vid)
		if pos == 'UpperRight':
			vid = int(vid) + 1
		if pos == 'LowerLeft':
			vid = int(vid)
	else:
		vid = int(vid)

	if x_res < CeilLen/2.0:
		print('xxxx')
		if pos == 'UpperLeft':
			hid = int(hid)
		if pos == 'LowerRight':
			hid = int(hid) - 1
		if pos == 'UpperRight':
			hid = int(hid) - 1	
		if pos == 'LowerLeft':
			hid = int(hid)
	else:
		hid = int(hid)	
# 	print(hid, vid)					
	strhid = str(hid)
	while len(strhid) < 2:
		strhid = '0' + strhid

		
	strvid = str(vid)
	while len(strvid) < 2:
		strvid = '0' + strvid
	
	tile = 'h' + strhid + 'v' + strvid

	return hid, vid, tile	

#-----------------------------------------------------------------------
def get_cord_Sinusoidal(cord, numCeil = 1200, debug = 0):
	'''
	Function to get the geographical coordinates of the given region. 
	The coordinates are correspoding to Sinusodial grid.
	
	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates 
			   (Top Latitude, Bottom Latitude, Left Longitude, Right Longitude)
		
		numCeil : number of pixel in one tile, default value 1200
		
	Return
    ----------
		meshLat : array
				  longitude 
		meshLon : array
				  longitude 
		hidMin	: int
    		      vetical index of the point 
		hidMax	: int
    		  	 vetical index of the point 
		vidMin  : int
    		  	  vetical index of the point 
		vidMax  : int
    		  	  vetical index of the point
	'''
	import numpy as np
	
	UpperLeft = (cord[0], cord[2])
	UpperRight = (cord[0], cord[3])
	LowerLeft = (cord[1], cord[2])
	LowerRight = (cord[1], cord[3])

	tileInfor = []

	tileInfor.append(get_point_tile_Sinusoidal((UpperLeft), 'UpperLeft'))
	tileInfor.append(get_point_tile_Sinusoidal((UpperRight), 'UpperRight'))
	tileInfor.append(get_point_tile_Sinusoidal((LowerRight), 'LowerRight'))
	tileInfor.append(get_point_tile_Sinusoidal((LowerLeft), 'LowerLeft'))

	hid = []
	vid = []
	for item in tileInfor:
		hid.append(item[0])
		vid.append(item[1])
	
	hidMax = np.max(hid)
	hidMin = np.min(hid)

	vidMax = np.max(vid)
	vidMin = np.min(vid)

	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1

	# update the hid and vid
	hid = np.arange(hidMin, hidMax + 1, 1)
	vid = np.arange(vidMin, vidMax + 1, 1)
	
	GridDim = (num_v * numCeil, num_h * numCeil)
	meshLat = np.full(GridDim, np.nan)
	meshLon = np.full(GridDim, np.nan)
	
	if debug == 1:
		print(hid)
		print(vid)
	
	for hh in hid:
		for vv in vid:
			
			strhid = str(hh)
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(vv)
			while len(strvid) < 2:
				strvid = '0' + strvid
	
			tile = 'h' + strhid + 'v' + strvid
			
			print(tile)
			if debug == 1:
				print('\n',tile)
		
			latitude, longitude = cal_sinu_grid(tile, numCeil)
			if debug == 1:
				print(np.max(latitude), np.min(latitude))
			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			
			if debug == 1:
				print(hIdx, vIdx)
				print(vIdx * numCeil, (vIdx + 1) * numCeil, hIdx * numCeil, (hIdx + 1) * numCeil)
		
			meshLat[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = latitude
				
			meshLon[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = longitude	
	return meshLat, meshLon, hidMin, hidMax, vidMin, vidMax

#-----------------------------------------------------------------------
'''
Function of PlateCarree...
'''
def get_point_tile_PlateCarree(cord, pos):
	
	lat = cord[0]
	lon = cord[1]

# 	print('  - get_point_tile_VNP46A1', lat // 10, lon // 10)
	vid = 8 - lat // 10
	hid = 18 + lon // 10
	
	res_lat = lat%10
	res_lon = lon%10
	
	if res_lat == 0:
		if pos == 'UpperLeft':
			vid = int(vid) + 1
		if pos == 'LowerRight':
			vid = int(vid)
		if pos == 'UpperRight':
			vid = int(vid) + 1
		if pos == 'LowerLeft':
			vid = int(vid)

	if res_lon == 0:
		if pos == 'UpperLeft':
			hid = int(hid)
		if pos == 'LowerRight':
			hid = int(hid) - 1
		if pos == 'UpperRight':
			hid = int(hid) - 1	
		if pos == 'LowerLeft':
			hid = int(hid)
	
	strhid = str(hid)
	while len(strhid) < 2:
		strhid = '0' + strhid

	strvid = str(vid)
	while len(strvid) < 2:
		strvid = '0' + strvid
		
	tile = 'h' + strhid + 'v' + strvid

	return hid, vid, tile

#-----------------------------------------------------------------------
def cal_PlateCarree_grid(tile, numCeil):
	
	import numpy as np
	
	vid = np.int(np.float(tile[4:6]))
	hid = np.int(np.float(tile[1:3]))
	print('  - Calulating geographical coordinates of ', tile)
	print('    Vertical Tile:', vid, 'Horizontal Tile:', hid)	
	
	latBoundary = [(8 - vid) * 10, (9 - vid) * 10]
	lonBoundary = [(hid - 18) * 10, (hid - 17) * 10]
	
	latitude  = np.linspace(latBoundary[1], latBoundary[0], numCeil)

	longitude  = np.linspace(lonBoundary[0], lonBoundary[1], numCeil)

	latitude = (latitude * np.ones((numCeil,1),np.float32)).T
	
	longitude = np.ones((numCeil,1),np.float32) * longitude
    
	return latitude, longitude

#-----------------------------------------------------------------------
def get_cord_PlateCarree(cord, numCeil = 2400, debug = 0):
	'''
	Function to get the geographical coordinates of the given region. 
	The coordinates are correspoding to Sinusodial grid.
	
	Parameters
    ----------
    	cord : list or tuple
			   geographical coordinates 
			   (Top Latitude, Bottom Latitude, Left Longitude, Right Longitude)
		
		numCeil : number of pixel in one tile, default value 1200
		
	Return
    ----------
		meshLat : array
				  longitude 
		meshLon : array
				  longitude 
		hidMin	: int
    		      vetical index of the point 
		hidMax	: int
    		  	 vetical index of the point 
		vidMin  : int
    		  	  vetical index of the point 
		vidMax  : int
    		  	  vetical index of the point
	'''
	import numpy as np

	UpperLeft = (cord[0], cord[2])
	UpperRight = (cord[0], cord[3])
	LowerLeft = (cord[1], cord[2])
	LowerRight = (cord[1], cord[3])

	tileInfor = []

		
	tileInfor.append(get_point_tile_PlateCarree(UpperLeft, 'UpperLeft'))
	tileInfor.append(get_point_tile_PlateCarree(UpperRight, 'UpperRight'))
	tileInfor.append(get_point_tile_PlateCarree(LowerRight, 'LowerRight'))
	tileInfor.append(get_point_tile_PlateCarree(LowerLeft, 'LowerLeft'))

	hid = []
	vid = []
	for item in tileInfor:
		hid.append(item[0])
		vid.append(item[1])
	
	hidMax = np.max(hid)
	hidMin = np.min(hid)

	vidMax = np.max(vid)
	vidMin = np.min(vid)
	
	
	num_h = hidMax - hidMin + 1
	num_v = vidMax - vidMin + 1

	# update the hid and vid
	hid = np.arange(hidMin, hidMax + 1, 1)
	vid = np.arange(vidMin, vidMax + 1, 1)
	
	GridDim = (num_v * numCeil, num_h * numCeil)
	meshLat = np.full(GridDim, np.nan)
	meshLon = np.full(GridDim, np.nan)
	
	if debug == 1:
		print('\n  - get_cord_VNP46A1 - hid: ', hid)
		print('\n  - get_cord_VNP46A1 - vid: ', vid)
	

	for hh in hid:
		for vv in vid:	
			strhid = str(hh)
			while len(strhid) < 2:
			
				strhid = '0' + strhid

			strvid = str(vv)
			while len(strvid) < 2:
				strvid = '0' + strvid
	
			tile = 'h' + strhid + 'v' + strvid
	
			if debug == 1:
				print('\n  - get_cord_VNP46A1 - tile: ',tile, hh, vv)
			
			latitude, longitude = cal_PlateCarree_grid(tile, numCeil)
			
			
			if debug == 1:
				print('\n  - get_cord_VNP46A1 - latitude: ',latitude[0,0], latitude[-1,0])
			hIdx = hh - hidMin 
			vIdx = vv - vidMin
			
			if debug == 1:
				print(hIdx, vIdx)
				print(vIdx * numCeil, (vIdx + 1) * numCeil, hIdx * numCeil, (hIdx + 1) * numCeil)
		
			meshLat[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = latitude
				
			meshLon[vIdx * numCeil : (vIdx + 1) * numCeil, \
					hIdx * numCeil : (hIdx + 1) * numCeil, ] = longitude	
					
	return meshLat, meshLon, hidMin, hidMax, vidMin, vidMax
	
#-----------------------------------------------------------------------
def plateCarree_to_geog(cord):
	'''
	Function of converting the Plate Carree to geographical projection
	https://en.wikipedia.org/wiki/Equirectangular_projection
	
	Parameters
    ----------
		sinusoidal point - list or tuple liked, (x, y)
		
	Return
    ----------
		geographical coordinates - latitude, longituede
	
	
	'''
	import numpy as np
	from numpy import sin, cos
	x = cord[0]
	y = cord[1]
	pi = 180.0 / np.pi
	R = 6371007.181000
	
	lamda_0 = 0.0
	phi_1 = 0.0
	
	lamda = x / R / cos(phi_1) + lamda_0
	phi	  = y / R  + phi_1

	lat = phi * pi
	lon = lamda * pi

	return lat, lon

#-----------------------------------------------------------------------
def geog_to_plateCarree(cord):
	'''
	Function of converting the geographical projection to Plate Carree
	https://en.wikipedia.org/wiki/Equirectangular_projection
	
	Parameters
    ----------
		geographical coordinates - list or tuple liked, (latitude, longituede)
		
	Return
    ----------
		equirectangular point - (x, y)
	
	
	'''
	import numpy as np
	from numpy import sin, cos
	
	pi = 180.0 / np.pi
	R = 6371007.181000
	
	lamda_0 = 0.0
	phi_1 = 0.0
	
	lat = cord[0]
	lon = cord[1]
	
	phi = lat / pi
	lamda = lon / pi	
	
	x = R * (lamda - lamda_0) * cos(phi_1)
	y = R * (phi - phi_1)
	
	return x, y	

#-----------------------------------------------------------------------
def get_tile_boundary_plateCarree(tile):

	import numpy as np
	vid = np.int(np.float(tile[4:6]))
	hid = np.int(np.float(tile[1:3]))
	print(' - Calulating geographical coordinates of ', tile)
	print(' - Vertical Tile:', vid, 'Horizontal Tile:', hid)	
	
	
	south  = (8 - vid) * 10.0
	north = (9 - vid) * 10.0
	west   = (hid - 18) * 10
	east   = (hid - 17) * 10

	return north, south, west, east

#-----------------------------------------------------------------------
def get_tiles(cord):

	import numpy as np
	
	hid_top, vid_top, _ = get_point_tile_PlateCarree((cord[0], cord[2]), pos = 'UpperLeft')
	hid_bot, vid_bot, _ = get_point_tile_PlateCarree((cord[1], cord[3]), pos = 'LowerRight')

# 	print( hid_top, vid_top )
# 	print( hid_bot, vid_bot ) 

	hids = np.arange(hid_top, hid_bot + 1, 1)
	vids = np.arange(vid_top, vid_bot + 1, 1)

	tiles = []
	for hid in hids:
		for vid in vids:
			strhid = str(hid)
			while len(strhid) < 2:
				strhid = '0' + strhid

			strvid = str(vid)
			while len(strvid) < 2:
				strvid = '0' + strvid
		
			tile = 'h' + strhid + 'v' + strvid
			tiles.append(tile)
	return tiles
	
#-----------------------------------------------------------------------	
def get_city_cords(cityname):

	cords = {}
	cords['Los Angeles']     = [35,   33,   -119,   -117]
	cords['Chicago']         = [42.5, 41.5, -88.5,  -87]
	cords['Denver']          = [40.5, 39.5, -105.5, -104.5]
	cords['San Fracisco']    = []
	cords['New York']        = [42, 40, -73, -71]
	cords['Seattle']         = [48, 47, -123, -122]
	cords['Washington D.C.'] = []
	cords['Boston']          = []
	cords['Houston']         = [30.25, 29.25,  -96,  -94.8]
	
	return cords[cityname]



#-----------------------------------------------------------------------	
def alber_equal_area(lon, lat, lat_0 = 40, lat_1 = 20, lat_2 = 60, lon_0 = -96):
	'''
	Function to convert ... to alber equal area
	'''
	
	import numpy as np
	from numpy import cos, sin, log

	def cal_alpha(phi, e):

		def first_term(phi, e):
	
			tan_phi = sin(phi) / ( 1 - (e*sin(phi))**2 )
	
			return tan_phi

		def second_term(phi, e):

			xx = ( 1 / (2 * e) ) * log( (1 - e * sin(phi) ) / ( 1 + e * sin(phi) ) )
	
			return xx

		alpha = (1 - e**2) * ( first_term(phi, e) - second_term(phi, e) )

		return alpha
		
	pi = 180.0 / np.pi
	R = 6378137
	# flattening
	f = 1/298.257233
	# eccentricity
	e = (2*f - f**2)**0.5
	
	lamda = lon / pi
	
	phi = lat / pi

	phi_0 = lat_0 / pi
	
	phi_1 = lat_1 / pi
	
	phi_2 = lat_2 / pi
	
	lamda_0 = lon_0 / pi
	
	alpha = cal_alpha(phi, e)
	
	alpha_0 = cal_alpha(phi_0, e)
	
	alpha_1 = cal_alpha(phi_1, e)
	
	alpha_2 = cal_alpha(phi_2, e)
	
	m1 = cos(phi_1) / ( 1 - (e * sin(phi_1))**2)**0.5
	
	m2 = cos(phi_2) / ( 1 - (e * sin(phi_2))**2)**0.5
	
	n = (m1**2 - m2**2) / (alpha_2 - alpha_1)
	
	C = m1**2 + n * alpha_1
	
	theta = n * (lamda - lamda_0)

	rho = R/n * (C - n*alpha)**0.5
	
	rho_0 =  R/n * (C - n*alpha_0)**0.5

	x = rho * sin(theta)
	
	y = rho_0 - rho*cos(theta)
	
	return x, y
	
