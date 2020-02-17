#----------------------------------------------------------------------------
def getWeights(cPixel, surPixels, obs):
	'''
	cPixel,    [dictionary-like], contains latitude, longitude of a standard 
			   pixel  
	surPixels, [dictionary-like], contains latitude, longitude of a raw pixel 
			   surround the center pixel
	obs,       [array-like], data to be resampled

	Albers Equal Area Projection is used in this function to calculate
	the overlap aera of the pixel of interest and the area

	'''
	import numpy as np
	from shapely.geometry import Polygon
	
	o = 1e-3
	cCord = [
			 (cPixel['Upper_Left_Lon'], cPixel['Upper_Left_Lat']),\
			 (cPixel['Upper_Right_Lon'], cPixel['Upper_Right_Lat']),\
		     (cPixel['Lower_Right_Lon'], cPixel['Lower_Right_Lat']),\
		     (cPixel['Lower_Left_Lon'], cPixel['Lower_Left_Lat'])
		    ]	
	cPolygon = Polygon(cCord)
	
# 	print('  - total area: %8.3f: '%(cPolygon.area))
	
	area = np.full_like(surPixels['Upper_Left_Lon'], np.nan)

	surPolygons = []
	
	Upper_Left_Lon	= surPixels['Upper_Left_Lon']
	Upper_Left_Lat	= surPixels['Upper_Left_Lat']
	Upper_Right_Lon = surPixels['Upper_Right_Lon']
	Upper_Right_Lat = surPixels['Upper_Right_Lat']
	Lower_Right_Lon = surPixels['Lower_Right_Lon']
	Lower_Right_Lat = surPixels['Lower_Right_Lat']
	Lower_Left_Lon 	= surPixels['Lower_Left_Lon']
	Lower_Left_Lat 	= surPixels['Lower_Left_Lat']
	
	
	for i in range(len(area[:,0])):
		for j in range(len(area[0,:])):
			#  modified here
			Cord = [(Upper_Left_Lon[i,j], Upper_Left_Lat[i,j]),\
					(Upper_Right_Lon[i,j], Upper_Right_Lat[i,j]),\
					(Lower_Right_Lon[i,j], Lower_Right_Lat[i,j]),\
					(Lower_Left_Lon[i,j], Lower_Left_Lat[i,j])]
			surPolygon = Polygon(Cord)
	    				 
			area[i,j] = surPolygon.intersection(cPolygon).area
			surPolygons.append(surPolygon)
			
	cArea = cPolygon.area
	overlapArea = np.nansum(area)
	if overlapArea < o:
		weights = np.full_like(surPixels['Upper_Left_Lon'], np.nan)
	else:
		weights = np.true_divide(area, np.nansum(area))


	return [weights, area, cArea]

#----------------------------------------------------------------------------	
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
	
#----------------------------------------------------------------------------	
