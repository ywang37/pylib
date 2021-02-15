#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Meng Zhou

Package of plotting satellite observation into 2D image...

To avoid library confliction, all libraries are privately import...

"""
#-----------------------------------------------------------------------
def define_map(cord, **kwargs):
	'''
	Function to define a geo-reference figure
	
	------------
	Parameters:
	cord : list, [North, South, West, East]
	
	------------
	Returns:
	fig	 : figure instance 
	ax	 : axes instance  
	proj : projection method
	
	'''
	
	import matplotlib.pyplot as plt
	import cartopy.crs as ccrs
	import numpy as np
	# total figure size
	figsize = kwargs.get('figsize', (12, 8))	
	# actual canvas that hold the figure
	axsize  = kwargs.get('axsize', [0.15,0.15,0.7,0.7])
	# define the central of the map...
	central_longitude  = kwargs.get('central_longitude', 0)
	# projection method
	proj    = kwargs.get('proj', ccrs.PlateCarree(central_longitude=central_longitude))
	
	base_proj = ccrs.PlateCarree(central_longitude=central_longitude)
	
	topLat = cord[0]
	botLat = cord[1]
	lefLon = cord[2]
	rigLon = cord[3]
	

	fig = plt.figure(figsize=figsize)

	ax = fig.add_axes(axsize, projection=proj)

	ax.set_extent([cord[2], cord[3], cord[0], cord[1]], crs=base_proj)
# 	ax.set_xlim(cord[2], cord[3])
# 	ax.set_ylim(cord[1], cord[0])	
	return [fig, ax, proj]
	

#-----------------------------------------------------------------------
def multiFigure(nRow, nCol, **kwargs):
	'''
	This function creates multi-figure layout...
	
	Parameter:
	nRow	: row number of the figure
	nCol	: column number of the figure
	
	Optional Paras:
	proj	: projection method (current only allow one projection method
	nUnit	: number of grid per figure
	figsize : figure size
	nGap	: number of grid between each figure
	nPlot	: number of plot want to display on the figure
	
	Output:
	
	figure  : figure instance 
	axes	: list of axes instance
	proj	: projection method (current only allow one projection method)
	
	here I suggest to use 9 * 9 for single figure, the corresponding nUint
	should be set as 10, for multiple figures, keep the default nUint
	
	'''
	import cartopy.crs as ccrs
	import matplotlib.pyplot as plt
	import matplotlib.gridspec as gridspec
	import matplotlib.patheffects as path_effects
	
	proj = kwargs.get('proj', ccrs.PlateCarree())
	nUint = kwargs.get('nUint', 25)
	nColGap = kwargs.get('nColGap', 5)
	nRowGap = kwargs.get('nRowGap', 5)
	nPlot = kwargs.get('nPlot', None)

	figsize = kwargs.get('figsize', (9, 9))
	projPos = kwargs.get('projPos', [])
	fontsize = kwargs.get('fontsize', 22)
	xlabel = kwargs.get('xlabel', 0.1)
	ylabel = kwargs.get('ylabel', 0.975)
	numOn = kwargs.get('numOn', True)
	
	cord = kwargs.get('cord', [90, -90, -180, 180])
	
	nRow_grid = nRow * nUint #+ nRowGap * (nCol - 1)
	nCol_grid = nCol * nUint #+ nColGap * (nCol - 1)
	
	gs = gridspec.GridSpec(nRow_grid,nCol_grid)

	axes = []
	fig = plt.figure(figsize = figsize)

	if nPlot == None:
		nPlot = nRow * nCol - 1
	else:
		nPlot = nPlot - 1
	
	for i in range(nRow):
		for j in range(nCol):
			nFigure = i * nCol + j
			if nFigure <= nPlot:
				txt_number = '(' + chr(97 + nFigure) + ')'		
				if nFigure in projPos:
# 					print(i*nUint, (i + 1)*nUint - nGap, j*nUint, (j + 1)*nUint - nGap)
# 					print(i, (i + 1)*nUint, nColGap, j, (j + 1)*nUint, nRowGap)
# 					print(i*nUint, (i + 1)*nUint - nColGap, j*nUint, (j + 1)*nUint - nRowGap)
# 					print()
					instant = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nColGap,  \
												 j*nUint:(j + 1)*nUint - nRowGap], \
												 projection=proj)
					instant.set_extent([cord[2], cord[3], cord[0], cord[1]])
					if numOn:            
						text  = instant.annotate(txt_number, xy=(xlabel, ylabel), 
												 xytext = (xlabel, ylabel), \
												 textcoords='axes fraction', \
												 xycoords=proj._as_mpl_transform(instant), \
												 color = 'k', ha='right', va='top', \
												 fontsize = fontsize)
						text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w'), \
											   path_effects.Normal()])
				else:
					instant = fig.add_subplot(gs[i*nUint:(i + 1)*nUint - nColGap,  \
												 j*nUint:(j + 1)*nUint - nRowGap]  )
					if numOn:
						text = instant.annotate(txt_number, xy=(xlabel, ylabel), \
												xytext = (xlabel, ylabel), \
												textcoords='axes fraction', \
												color = 'k', ha='right', va='top',\
												fontsize = fontsize)								 
						text.set_path_effects([path_effects.Stroke(linewidth=2, foreground='w'), \
											   path_effects.Normal()])
				axes.append(instant)

	return [fig, axes, proj]

#-----------------------------------------------------------------------
def save_figure(fig, saveName, dpi = 300):
	'''
	Function to save a figure object...
	
	'''

	import matplotlib.pyplot as plt
	
	print(' - Saving', saveName)
	
	fig.savefig(saveName, bbox_inches='tight', dpi=dpi)
	plt.close()

#-----------------------------------------------------------------------
def plotSat_2d_basic(fig, ax, proj, plotData, channel, **kwargs):

	'''
	Function to create the basic basemap projection figure...
    
	Parameters
    ----------    
	fig	 	 : figure instance 
	ax	 	 : axes instance  
    proj     : projection method
    plotData : list of dictionaries, it holds the data to be plotted on the map
    		   fixed-key: 'latitude', 'longitude', array, hold the coordinates
			   user-defined key: array, holds the observational data...   
    channel	 : key of the data to be plotted that are named in plotData dictionary
        
	------------
	Parameters:
	fig	 	 : figure instance 
	ax	 	 : axes instance  
	cbar_pos : orientation of color bar, vertical or horizontal
	
	Other optional parameters are explained in the code
	
	------------
	Returns:
	ax	 	 : axes instance 
	
	'''

	import numpy as np
	import matplotlib.pyplot as plt
	import cartopy.crs as ccrs
	import cartopy.feature as cfeature
	import matplotlib.ticker as ticker
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
	

	def cord_Compare(currentMax, currentMin, data):
	
		if currentMax < np.max(data):
			currentMax = np.max(data)
		
		if currentMin > np.min(data):
			currentMin = np.min(data)
	
		return currentMax, currentMin
	
	# default parameters are set here
	# range of the value [min, max]
	vRange		= kwargs.get('vRange', [0,1])
	# title of the figure
	title		= kwargs.get('title', '')
	
	# coordinate of the figure, [North, South, West, East], if it is not
	# provided, the code will use the private function cord_Compare
	# generate the cord
	cord		= kwargs.get('cord', False)
	
	scaleFactor	= kwargs.get('scaleFactor', 1.0)
	
	# line width of lines in the figure, including river, costaline etc.
	linewidth  	= kwargs.get('linewidth',1)
	fontsize  	= kwargs.get('fontsize', 12)
	lineColor	= kwargs.get('lineColor', 'dimgray')
	
	# relate to the color of the figure...
	# color map of the figure, string
	cmap		= kwargs.get('cmap', 'jet')
	# numer of the color used in the figure, 255 denotes a continuous
	# color bar
	numClr   	= kwargs.get('numClr', 255)
	
	# specify bad value of the data, if None, it will just mask out the 
	# NaN value
	badValue	= kwargs.get('badValue', None)
	# color used to denote the bad value
	badClr 	    = kwargs.get('badClr', None)
	
	# color used to denote pixel whose value is smaller than minimum	
	underClr	= kwargs.get('underClr', None)
	# color used to denote pixel whose value is greater than max	
	overClr	    = kwargs.get('overClr', None)
	# color of the background
	bgClr       = kwargs.get('bgClr', None)
	# color of waterbody
	waterClr    = kwargs.get('waterClr', None)
	
	
	# relate to axis ticks...
	# specify the direction of the label, refer to
	# https://scitools.org.uk/cartopy/docs/latest/gallery/tick_labels.html
	dateline_direction_label = kwargs.get('dateline_direction_label', False)
	# specify the symbol of degree, default one is upper circle 
	degree_symbol = kwargs.get('degree_symbol', u'\u00B0')
	# format of the coordinates
	number_format        = kwargs.get('number_format', '4.2f') 
	# specify the direction of the label, refer to
	# https://scitools.org.uk/cartopy/docs/latest/gallery/tick_labels.html
	zero_direction_label = kwargs.get('zero_direction_label', False)
	# number of ticks of the latitude and longitude 
	numTicks_Y           = kwargs.get('numTicks_Y', 5)
	numTicks_X           = kwargs.get('numTicks_X', 5)
	# logical switch of show the latitude or longitude ticks
	yvisible 	         = kwargs.get('yvisible', True)
	xvisible 	         = kwargs.get('xvisible', True)
	# if not None, it will be used as the label of the X and Y tick
	ticks_X  	         = kwargs.get('ticks_X',  None)
	ticks_Y  	         = kwargs.get('ticks_Y',  None)
	
	# related to color bar..
	cbarOn      = kwargs.get('cbarOn', True)
	# orientation of color bar, vertical or horizontal
	cbar_pos    = kwargs.get('orientation', 'vertical') # vertical or horizontal
	
	
	data_proj = ccrs.PlateCarree()
	
	#--------------
	# initialize the maximum and minimum of the latitude and longitude
	maxLat = 0
	minLat = 0
	maxLon = 0
	minLon = 0
	
	# get the color range (maximum, minimum) of the data
	vmin = vRange[0]
	vmax = vRange[1]
	
	#  - set the color map
	cmap = plt.get_cmap(cmap, numClr)
	
	if underClr is not None:
		cmap.set_under(underClr)

	if overClr is not None:
		cmap.set_over(overClr)

	# set the color for the NaN value
	if badClr is not None:
		cmap.set_bad(color = badClr)
	
	if waterClr is None:
		waterClr = bgClr
	
	# we plot the figure by orbit...
	for i in range( len(plotData)):

		lat = plotData[i]['latitude']
		lon = plotData[i]['longitude']

		lat = calPixelEdge(lat)
		lon = calPixelEdge(lon)
		
		plotChannel = plotData[i][channel][:] * scaleFactor

		if badValue is None:
			plotChannel = np.ma.masked_array(plotChannel, plotChannel != plotChannel)			
		else:
			plotChannel[np.where(plotChannel == badValue)] = np.nan
			plotChannel = np.ma.masked_array(plotChannel, plotChannel != plotChannel)			


		maxLat, minLat = cord_Compare(maxLat, minLat, lat)
		maxLon, minLon = cord_Compare(maxLon, minLon, lon)
		#----------------------------
		# data is actually plot here
		#----------------------------

		sc = ax.pcolormesh(lon, lat, plotChannel, transform = data_proj, \
		                   cmap = cmap, vmin = vmin, vmax = vmax)
		                     
	# plot the buildin coastline, country, state lines
	# create the label of the x and y axis
	if	cord == False:
		labelLat = np.linspace(minLat, maxLat, numTicks_Y)
		labelLon = np.linspace(minLon, maxLon, numTicks_X)
		minor_labelLat = np.linspace(minLat, maxLat, 4*numTicks_Y)
		minor_labelLon = np.linspace(minLon, maxLon, 4*numTicks_X)
	else:
		labelLat = np.linspace(cord[0], cord[1], numTicks_Y)
		labelLon = np.linspace(cord[2], cord[3], numTicks_X)
		minor_labelLat = np.linspace(cord[0], cord[1], 4*numTicks_Y)
		minor_labelLon = np.linspace(cord[2], cord[3], 4*numTicks_X)

	
	if ticks_X is not None:
		labelLon = ticks_X
	if ticks_Y is not None:
		labelLat = ticks_Y
	
	
	# set the paralle line of the latitude and longitude	
	if xvisible == True:
		ax.set_xticks(labelLon, crs=data_proj)
# 		ax.set_xticks(minor_labelLon, minor = True, crs=data_proj)
		lon_formatter = LongitudeFormatter(
				zero_direction_label=zero_direction_label,
				dateline_direction_label=dateline_direction_label,
				number_format=number_format, degree_symbol=degree_symbol)
		ax.xaxis.set_major_formatter(lon_formatter)
		
	if yvisible == True:
		ax.set_yticks(labelLat, crs=data_proj)	
# 		ax.set_yticks(minor_labelLat, minor = True, crs=data_proj)		
		lat_formatter = LatitudeFormatter(number_format=number_format,
				degree_symbol=degree_symbol)		
		ax.yaxis.set_major_formatter(lat_formatter)

	
	ax.gridlines(linewidth=linewidth*0.75, crs=data_proj, xlocs=labelLon, ylocs=labelLat, 
	             edgecolor = lineColor,linestyle=':')
	
	# set the fontsize of the X & Y axis...
	ax.tick_params( labelsize = fontsize)
	
	# add the state boundary here
	ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor)
	# add the lake boundary here
	ax.add_feature(cfeature.LAKES.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor, alpha=0.5, color = waterClr)
	# add the river line here
	ax.add_feature(cfeature.RIVERS.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor)
	               
	# add the coastline, including major islands.
	ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor)
	# add the land polygons, including major islands.
	ax.add_feature(cfeature.LAND.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor,color = bgClr)
	# add the country boundaries.
	ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=linewidth, \
	               edgecolor = lineColor,linestyle=':')
	# add the ocean
	ax.add_feature(cfeature.OCEAN, color = waterClr)

	# also set the background to the background color	
	ax.background_patch.set_facecolor(bgClr)
	
	# then set the color bar...
	if cbarOn == True:
		# set the color bar
		if cbar_pos == "horizontal":		
			# x_off = 0, y_off= -0.1 are suggested here
			set_colorbar(fig, ax, cbar_pos, **kwargs)
		if cbar_pos == "vertical":	
			# x_off = 0.015,  y_off= 0 are suggested here
			set_colorbar(fig, ax, cbar_pos, **kwargs)
			
	# set the title
	ax.set_title(title, fontsize = fontsize)
	return ax
	
#-----------------------------------------------------------------------
def plotSat_2d_scatter(ax, data, **kwargs):
	'''
	'''
	import numpy as np
	import matplotlib as mpl
	import matplotlib.pyplot as plt
	
	# default parameters are set here
	vRange		= kwargs.get('vRange', [0,1])
	lineColor	= kwargs.get('lineColor', 'w')
	cmap		= kwargs.get('cmap', 'jet')
	cmapSTD		= kwargs.get('cmapSTD', 'Wistia')
	badColor	= kwargs.get('badColor', 'gray')
	marker 		= kwargs.get('marker', ".")
	markersize	= kwargs.get('markersize', 200)
	markersizeSTD	= kwargs.get('markersizeSTD', 400)
	color		= kwargs.get('color', 'red')
	edgecolors	= kwargs.get('edgecolors', None)
	linewidths	= kwargs.get('linewidths', None)
	colorMap	= kwargs.get('colorMap', False)
	plotValue	= kwargs.get('plotValue', False)
	plotSTD		= kwargs.get('plotSTD', False)
	fmt			= kwargs.get('fmt', 'o')
	elinewidth	= kwargs.get('elinewidth', 1)
	capsize		= kwargs.get('capsize', 2)
	numColr		= kwargs.get('numColr', 255)
	absentColor = kwargs.get('absentColor', None)
	
	yerrs		= kwargs.get('yerr', 0)
	yerrs		= kwargs.get('xerr', 0)
	scale		= kwargs.get('xerr', 20)
		
	xs = np.array(data['longitude'])
	ys = np.array(data['latitude'])

	# set the color map
	cmap = plt.get_cmap(cmap, numColr)
	
	if absentColor is not None:
		cmap.set_under(absentColor)

	
	if plotValue:
		zs = np.array(data['observation'])
		idx = np.where( zs == zs )
		zs = zs[idx]
		xs = xs[idx]
		ys = ys[idx]
		if plotSTD:
			std = np.array(data['std'])
			xerrs = std[idx]
		
	vmin = vRange[0]
	vmax = vRange[1]
	
	t = mpl.markers.MarkerStyle(marker=marker)
	t._transform = t.get_transform().rotate_deg(45)
		
	
	if (plotValue == True) & (plotSTD == False) :
		sc = ax.scatter(xs, ys, s = markersize, c = zs, marker=t, cmap=cmap, vmin = vmin,\
						vmax = vmax, linewidths = linewidths, edgecolors = edgecolors)
						
	if (plotValue == True) & (plotSTD == True):
		ax.scatter(xs, ys, s = markersizeSTD, c = xerrs, marker=t, cmap=cmapSTD, vmin = 0,\
						vmax = 0.2, linewidths = linewidths, edgecolors = edgecolors)
								
		sc = ax.scatter(xs, ys, s = markersize, c = zs, marker=t, cmap=cmap, vmin = vmin,\
						vmax = vmax, linewidths = linewidths, edgecolors = edgecolors)


	if (plotValue == False) & (plotSTD == False):
		sc = ax.scatter(xs, ys, s = markersize, c = color, marker=t, linewidths = linewidths,\
		                edgecolors = edgecolors)
		                
	if colorMap == True:
		return ax, sc
	else:
		return ax

#-----------------------------------------------------------------------
def plotSat_2d_Traj(ax, lat, lon, **kwargs):

# 	color = kwargs.get('color', 'r')
	ax.plot(lon, lat, **kwargs)
	
	return ax

#-----------------------------------------------------------------------
def plotSat_1d_Prof(ax, data, channel, **kwargs):
	'''
	Function to plot the satellite profile...
	
	'''
	import numpy as np
	import matplotlib.pyplot as plt

	scaleFactor	= kwargs.get('scaleFactor', 1.0)
	vRange		= kwargs.get('vRange', None)
	cbar_label  = kwargs.get('cbar_label', '')
	title		= kwargs.get('title', '')
	y_label		= kwargs.get('y_label', '')
	x_label		= kwargs.get('x_label', '')
	cord		= kwargs.get('cord', False)
	lineWidth  	= kwargs.get('lineWidth', 1.5)
	fontsize  	= kwargs.get('fontsize', 12)
	lineColor	= kwargs.get('lineColor', 'dimgray')	
	stdOn	= kwargs.get('stdOn', False)	


	lat = data['latitude']
	lon = data['longitude']
	plotChannel = data[channel]
	plotChannel = np.ma.masked_array(plotChannel, plotChannel != plotChannel)
	
	xdata = np.zeros(len(plotChannel))
	for idx in range(len(plotChannel)):
		xdata[idx] = idx	
	
	if vRange == None:
		vRange = []
		vRange.append(np.min(plotChannel)-0.1)
		vRange.append(np.max(plotChannel)+0.1)
		
	if stdOn == True:
		error = data['std']
		
		ax.errorbar(xdata, plotChannel, yerr = error,  color = lineColor, ecolor = 'red', elinewidth = lineWidth/2)	
			
	ax.plot(xdata, plotChannel, color = lineColor, linewidth =lineWidth, label = channel)


	ax.set_title(title,fontsize=fontsize)
	ax.set_xlim(0, np.max(xdata)+.1)
	ax.set_ylim(vRange[0], vRange[1])
	ax.set_ylabel(y_label, fontsize=fontsize)

	numTick = 5
	lenX = len(plotChannel)
	idxIni = np.int((lenX - np.int(((lenX-1)/numTick)) * (numTick - 1))/2)
	ticks = np.arange(idxIni, lenX-1, (lenX - 1)/numTick)	
	
	ax.xaxis.set_ticks(ticks)
	xaxis_idx = [int(idx) for idx in np.arange(1, lenX-1, (lenX-1)/numTick )]
	xticklabels = [r"%.1f"%lat[i]+"\n" + r"%.1f"%lon[i] for i in xaxis_idx]
	ax.xaxis.set_ticklabels( xticklabels , fontsize=fontsize)
	
	x_off = - 7
	ax.annotate(r'Lat ($^\circ$N)',   xy=(1, 0), xytext=(20, x_off-fontsize), va='bottom',  
			xycoords='axes fraction', textcoords='offset points', fontsize = fontsize)
	ax.annotate(r'Lon ($^\circ$E)',   xy=(1, 0), xytext=(20, x_off-2*fontsize), va='bottom',  
			xycoords='axes fraction', textcoords='offset points', fontsize = fontsize)

	return ax
	
#-----------------------------------------------------------------------	
def plotSat_1d_Prof_rev(ax, data, channel, **kwargs):
	'''
	Function to plot the satellite profile...
	
	'''
	import numpy as np
	import matplotlib.pyplot as plt

	scaleFactor	= kwargs.get('scaleFactor', 1.0)
	vRange		= kwargs.get('vRange', None)
	yRange		= kwargs.get('yRange', None)
	cbar_label  = kwargs.get('cbar_label', '')
	title		= kwargs.get('title', '')
	y_label		= kwargs.get('y_label', '')
	x_label		= kwargs.get('x_label', '')
	cord		= kwargs.get('cord', False)
	lineWidth  	= kwargs.get('lineWidth', 1.5)
	fontsize  	= kwargs.get('fontsize', 12)
	lineColor	= kwargs.get('lineColor', 'dimgray')	
	stdOn	= kwargs.get('stdOn', False)	


	lat = data['latitude']
	lon = data['longitude']
	plotChannel = data[channel]
	plotChannel = np.ma.masked_array(plotChannel, plotChannel != plotChannel)
	
	ydata = np.zeros(len(plotChannel))
	for idx in range(len(plotChannel)):
		ydata[idx] = idx	
	
	if vRange == None:
		vRange = []
		vRange.append(np.min(plotChannel)-0.1)
		vRange.append(np.max(plotChannel)+0.1)
		
	if stdOn == True:
		error = data['std']
		
		ax.errorbar(xdata, plotChannel, yerr = error,  color = lineColor, ecolor = 'red', elinewidth = lineWidth/2)	
	
# 	print(plotChannel)
	ax.plot(plotChannel, lat, color = lineColor, linewidth =lineWidth, label = channel)
	
	ax.set_title(title,fontsize=fontsize-2)
	ax.set_ylim(yRange[0], yRange[1])
	ax.set_xlim(vRange[0], vRange[1])
	ax.set_xlabel(y_label, fontsize=fontsize)
	
	numTick = 5

	latTicks = np.linspace(yRange[0], yRange[1], numTick)
	
	lonTicks = []
	
	lat = lat[lat == lat]
	lon = lon[lon == lon]
	for latTick in latTicks:
		dis = np.abs(lat - latTick)
		
		idx = np.where(dis == np.min(dis))
		
		print(latTick,idx)
		lonTicks.append(lon[idx])
	
	lonTicks = np.array(lonTicks)

	print(lonTicks, latTicks)

	ax.yaxis.set_ticks(latTicks)
	yticklabels = [r"%.1f"%dislat+"\n" + r"%.1f"%dislon for dislat, dislon in zip(latTicks, lonTicks)]
	ax.yaxis.set_ticklabels( yticklabels , fontsize=fontsize)

	
	x_off = - 12
	ax.annotate(r'Lat ($^\circ$N)',   xy=(-0.2, 0), xytext=(0, x_off-fontsize), va='bottom',  
			xycoords='axes fraction', textcoords='offset points', fontsize = fontsize-2)
	ax.annotate(r'Lon ($^\circ$E)',   xy=(-0.2, 0), xytext=(0, x_off-2*fontsize), va='bottom',  
			xycoords='axes fraction', textcoords='offset points', fontsize = fontsize-2)

	return ax	
		
#-----------------------------------------------------------------------
def interpolation2D(data, winLen, debug = 0):	
	'''
	
	-----------
	Parameters:
	data  : 	array data to be interpolated
	winLen:	length of the window to implement the interplation
	
	-----------
	return:
	data  : array data after interpolation
	
	function used to fill the unvalid value in a array
	Naive average is used here, gaussian average can be also apply, however, need
	to be further investigated
	'''
	
	import inspect
	import numpy as np
	
	def retrieve_name(var):
			"""
			Gets the name of var. Does it from the out most frame inner-wards.
			:param var: variable to get name from.
			:return: string
			"""
			for fi in reversed(inspect.stack()):
				names = [var_name for var_name, var_val in fi.frame.f_locals.items() if var_val is var]
				if len(names) > 0:
					return names[0]
	
	def gauss2D(shape=(3,3),sigma=0.5):
		"""
		2D gaussian mask - should give the same result as MATLAB's
		fspecial('gaussian',[shape],[sigma])
		"""
		import numpy as np
		
		m , n = [(ss-1.)/2. for ss in shape]
		y,x = np.ogrid[-m:m+1,-n:n+1]
		h = np.exp( -(x*x + y*y) / (2.*sigma*sigma) )
		h[ h < np.finfo(h.dtype).eps*h.max() ] = 0
		sumh = h.sum()
		if sumh != 0:
			h /= sumh
		return h
		
	my_var_name = retrieve_name(data)
	print(' - Interpolating variable' , my_var_name, '\n')
	
	
	if winLen//2 == 0:
		winLen = winLen + 1
	
	invalid_idx = np.where(data != data)
	
	data[invalid_idx] = np.nan

	length, width = np.shape(data)
	halfWinLen = np.int((winLen - 1)/2)
	
	
	tempData = np.full((length + winLen - 1, width + winLen - 1), np.nan)
	tempData[halfWinLen:length + halfWinLen, halfWinLen: width + halfWinLen] = data
	
	tempData3D = np.full((length, width, winLen**2), np.nan)
	
	for i in range(winLen):
		for j in range(winLen):
			tempData3D[: , :, i*winLen + j] = tempData[i:i+length ,j:j+width]
            
	weights = gauss2D(shape=(winLen**2,1),sigma=1)
	weights = np.reshape(weights, len(weights))
	
	if debug == 1:
		print('weights', np.shape(weights))
		print('tempData3D', np.shape(tempData3D[0,0,:]))

# 	tempData2D = np.average(tempData3D, axis = 2, weights = weights)
	
	tempData2D = np.nanmean(tempData3D, axis = 2)
	
	if debug == 1:
		print(data[invalid_idx])
		print(tempData2D[invalid_idx])
	data[invalid_idx] = tempData2D[invalid_idx]
	
	return data
	
#-----------------------------------------------------------------------
def histeq(ims, nbr_bins=256 , debug = 0):
	
	'''
	im	: list of data...
	
	'''
	import numpy as np
	# get image histogram
	data = []
	for i in range(len(ims)):
		data.append(ims[i].flatten())
    	
	if debug == 1:
		print('Length of data: ',len(data))
    
	dataTemp = data[0]
	if len(data) > 1:
		for i in range(len(data)-1):
			dataTemp = np.concatenate((dataTemp, data[i+1]))
    			    
	imhist,bins = np.histogram(dataTemp[~np.isnan(dataTemp)],nbr_bins)
	cdf = imhist.cumsum()		# cumulative distribution function
	cdf = 255 * cdf / cdf[-1]	# normalize
    
    # use linear interpolation of cdf to find new pixel values
	im2 = []
	for i in range(len(data)):
		datatemp = np.interp(data[i], bins[:-1], cdf)
		im2.append( datatemp.reshape(ims[i].shape) )
		
	return im2, cdf
	
#-----------------------------------------------------------------------	
def histeq_multifigure(data_list, channels, nbr_bins=256,  debug = 0):
	
	import numpy as np
	print('\n - Calling histeq_multifigure')
	for channel in channels:
		print('\n - Do histogram equlization, Channel:', channel)
		data_total = []
		
		for data in data_list:
			print('\n - flatten', channel, 'data size',  len(data[channel].flatten()))
			data_total.append(data[channel].flatten())
		
		dataTemp = data_total[0]
		if len(data_total) > 1:
			for i in range(len(data_total)-1):
				dataTemp = np.concatenate((dataTemp, data_total[i+1]))		
				
		imhist,bins = np.histogram(dataTemp[~np.isnan(dataTemp)],nbr_bins)
		cdf = imhist.cumsum()		# cumulative distribution function
		cdf = 255 * cdf / cdf[-1]	# normalize

		for i in range(len(data_total)):
			datatemp = np.interp(data_total[i], bins[:-1], cdf)
			data_list[i][channel] = datatemp.reshape(data_list[i][channel].shape)
			
	return data_list

#-----------------------------------------------------------------------
def plot_truecolor(fig, ax, proj, plotData, channel, debug = 0, **kwargs):

	import numpy as np
	import matplotlib.pyplot as plt
	import cartopy.crs as ccrs
	import cartopy.feature as cfeature
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
	
	cord		= kwargs.get('cord', False)
	title		= kwargs.get('title', '')
	numTicks_Y  = kwargs.get('numTicks_Y', 5)
	numTicks_X  = kwargs.get('numTicks_X', 5)
	cord		= kwargs.get('cord', False)
	lineWidth  	= kwargs.get('lineWidth', 0.6)
	fontsize  	= kwargs.get('fontsize', 14)
	lineColor	= kwargs.get('lineColor', 'dimgray')
	numColr		= kwargs.get('numColr', 255)	
	absentColor = kwargs.get('absentColor', 'lightgray')
	dateline_direction_label = kwargs.get('dateline_direction_label', False) 
	degree_symbol = kwargs.get('degree_symbol', u'\u00B0') 
	number_format = kwargs.get('number_format', '4.2f') 
	zero_direction_label = kwargs.get('zero_direction_label', False)
	yvisible 	= kwargs.get('yvisible', True)
	xvisible 	= kwargs.get('xvisible', True)
	

	def cord_Compare(currentMax, currentMin, data):
	
		if currentMax < np.max(data):
			currentMax = np.max(data)
		
		if currentMin > np.min(data):
			currentMin = np.min(data)
	
		return currentMax, currentMin

	o = 1e-2

	maxLat = 0
	minLat = 0
	maxLon = 0
	minLon = 0
	
	for i in range(len(plotData)):

		lat = plotData[i]['latitude']
		lon = plotData[i]['longitude']
		plotChannel = plotData[i][channel]
		
		if debug == 1:
			print('  Debuging - Shape of lat: ', np.shape(lat))
			print('  Debuging - Shape of lon: ', np.shape(lon))
			print('  Debuging - Shape of plotChannel: ', np.shape(plotChannel))
										
		maxLat, minLat = cord_Compare(maxLat, minLat, lat)
		maxLon, minLon = cord_Compare(maxLon, minLon, lon)

		mesh_rgb = plotChannel[:, :, :]
		# seems here is a bug...
# 		mesh_rgb = plotChannel[:, :-1, :]

		colorTuple = mesh_rgb.reshape((mesh_rgb.shape[0] * mesh_rgb.shape[1]), 3)
		colorTuple = np.insert(colorTuple,3,1.0,axis=1)
			
		ax.pcolormesh(lon, lat, plotChannel[:,:,0], color = colorTuple,  transform=proj)
	

	# create the label of the x and y axis
	if	cord == False:
		labelLat = np.linspace(minLat, maxLat, numTicks_Y)
		labelLon = np.linspace(minLon, maxLon, numTicks_X)
	else:
		labelLat = np.linspace(cord[0], cord[1], numTicks_Y)
		labelLon = np.linspace(cord[2], cord[3], numTicks_X)
	
	# set the paralle line of the latitude and longitude, I personal prefer the line 
	# width is half width of other line...	
	
	if xvisible == True:

		ax.set_xticks(labelLon, crs=proj)
		lon_formatter = LongitudeFormatter(
				zero_direction_label=zero_direction_label,
				dateline_direction_label=dateline_direction_label,
				number_format=number_format, degree_symbol=degree_symbol)
		ax.xaxis.set_major_formatter(lon_formatter)
		
	if yvisible == True:
		ax.set_yticks(labelLat, crs=proj)
		lat_formatter = LatitudeFormatter(number_format=number_format,
				degree_symbol=degree_symbol)		
		ax.yaxis.set_major_formatter(lat_formatter)
		
	ax.tick_params( labelsize = fontsize)
	
	# add the states boundary here
	ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.5, edgecolor = lineColor)
		
	ax.background_patch.set_facecolor(absentColor)		
	
	# set the title of the figure
	ax.set_title(title, fontsize = fontsize)

	
	return ax

#-----------------------------------------------------------------------
def grouper(series, lenMin = 10):
	
	import numpy as np
	# - get the unique series
	unique_series = list(set(series))
	unique_series.sort()
	unique_series = np.array(unique_series)
	
	
	timeSeries = []
	
	for us in unique_series:
		timeSeries.append( (us//100) * 60 + (us%100))
	
	diff = np.diff(timeSeries)
	
	idx = np.array(np.where(diff>lenMin)) + 1
	
	idx =  np.concatenate((np.array([0]), np.reshape(idx, np.size(idx)), np.array([len(diff) + 2])))
	
	clutterSet = []
	
	for i in range(len(idx)-1):
		clutterSet.append(unique_series[idx[i]:idx[i+1]])
	
	re_set = {}
	
	for subSet in clutterSet:
		re_subset = []
		for i in range(len(subSet)):
			temp = str(subSet[i])
			
			while len(temp) < 4:
				temp = '0' + temp
			re_subset.append(temp)
			
		re_set[re_subset[0]] = re_subset
	
	return re_set
	
#-----------------------------------------------------------------------
def	plot_2d_gis(ax, proj, level, cord, **kwargs):

	'''
	Function to a high resolution GIS map as background...
	
	-----------
	Parameters:
	ax		: matplotlib axes instance
	proj	: Cartopy projection method
	level	: zoom-in level of a map, https://wiki.openstreetmap.org/wiki/Zoom_levels
	cord	: coordinate boundary of the figure, [North, South, West, East]
	-----------
	Returns:
	ax		: matplotlib axes instance
	
	
	Update by MZ, Cartopy 0.18 has fixed the issue of figure resolution,
	we don't need to specify the regrid_shape anymore after version 0.18 
	
	'''

	import numpy as np
	import matplotlib.pyplot as plt
	import cartopy.crs as ccrs
	import cartopy.feature as cfeature
	import cartopy.io.img_tiles as cimgt
	from owslib.wmts import _TILE_MATRIX_SET_TAG
	from owslib.wmts import _TILE_MATRIX_SET_LIMITS_TAG
	from owslib.wmts import _TILE_MATRIX_LIMITS_TAG
	from owslib.wmts import TileMatrixSetLink, TileMatrixLimits
	from time import sleep
	from random import randint
	
	# we define two functions here and overwrite the default function, since
	# the current default function may suffer some problem
	def custom_from_elements(link_elements):
		links = []
		for link_element in link_elements:
			matrix_set_elements = link_element.findall(_TILE_MATRIX_SET_TAG)
			if len(matrix_set_elements) == 0:
				raise ValueError('Missing TileMatrixSet in %s' % link_element)
			elif len(matrix_set_elements) > 1:
				set_limits_elements = link_element.findall(
					_TILE_MATRIX_SET_LIMITS_TAG)
				if set_limits_elements:
					raise ValueError('Multiple instances of TileMatrixSet'
									  ' plus TileMatrixSetLimits in %s' %
									  link_element)
				for matrix_set_element in matrix_set_elements:
					uri = matrix_set_element.text.strip()
					links.append(TileMatrixSetLink(uri))
			else:
				uri = matrix_set_elements[0].text.strip()

				tilematrixlimits = {}
				path = '%s/%s' % (_TILE_MATRIX_SET_LIMITS_TAG,
								  _TILE_MATRIX_LIMITS_TAG)
				for limits_element in link_element.findall(path):
					tml = TileMatrixLimits(limits_element)
					if tml.tilematrix:
						tilematrixlimits[tml.tilematrix] = tml

				links.append(TileMatrixSetLink(uri, tilematrixlimits))
		return links

	def new_get_image(self, tile):

		import six
		from PIL import Image
		if six.PY3:
			from urllib.request import urlopen, Request
		else:
			from urllib2 import urlopen, Request
	
		numSecond = np.float(randint(10,20))/10
		url = self._image_url(tile)
		print(' - Requesting ', url, 'waiting...', numSecond, 'seconds...')
		req = Request(url)
		sleep(numSecond)
		
		req.add_header('User-agent', 'your bot 0.1')
		fh = urlopen(req)

		im_data = six.BytesIO(fh.read())
		fh.close()
	
		img = Image.open(im_data)

		img = img.convert(self.desired_tile_form)
	# 	print(' - Function get_image:', np.shape(img))
	
		return img, self.tileextent(tile), 'lower'
	

	TileMatrixSetLink.from_elements = custom_from_elements
	cimgt.GoogleWTS.get_image = new_get_image

	title		= kwargs.get('title', '')
	linewidth  	= kwargs.get('linewidth', 1.5)
	fontsize  	= kwargs.get('fontsize', 14)
	linecolor	= kwargs.get('linecolor', 'dimgray')
	regrid_shape = kwargs.get('regrid_shape', (2500,2500))
	mstyle       = kwargs.get('mstyle', 'satellite')


	arcgis_url = 'https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}.jpg'
# 	arcgis_url ='https://mts0.google.com/vt/lyrs={style}@177000000&hl=en&src=api&x={x}&y={y}&z={z}&s=G'
	arcgis_url ='https://mts0.google.com/vt/lyrs={style}@177000000&hl=en&src=api&x={x}&y={y}&z={z}&s=G'
	tiles = cimgt.GoogleTiles(style=mstyle, url=arcgis_url)
	# 
# 	ax.add_image(tiles, level, regrid_shape = regrid_shape,\
# 	             interpolation='Gaussian')
	ax.add_image(tiles, level)
	ax = set_x_y_label(ax, proj, cord, **kwargs)
	
	# add the state line
# 	ax.add_feature(cfeature.STATES.with_scale('10m'), \
# 	               linewidth=linewidth, edgecolor = linecolor)
	
	# set title
	ax.set_title(title, fontsize = fontsize)
	
	return ax

#-----------------------------------------------------------------------
def set_x_y_label(ax, proj, cord, **kwargs):

	import numpy as np
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
	
	yvisible = kwargs.get('yvisible', True)
	xvisible = kwargs.get('xvisible', True)

	y_label		= kwargs.get('y_label', None)
	x_label		= kwargs.get('x_label', None)	
	fontsize  	= kwargs.get('fontsize', 14)
	ticks_X  	  = kwargs.get('ticks_X',  None)
	ticks_Y  	  = kwargs.get('ticks_Y',  None)

	
	dateline_direction_label = kwargs.get('dateline_direction_label', False) 
	degree_symbol = kwargs.get('degree_symbol', u'\u00B0') 
	number_format = kwargs.get('number_format', '4.2f') 
	zero_direction_label = kwargs.get('zero_direction_label', False)
	numTicks_Y  = kwargs.get('numTicks_Y', 5)
	numTicks_X  = kwargs.get('numTicks_X', 5)

	if ticks_X is not None:
		labelLat = ticks_X
	if ticks_Y is not None:
		labelLon = ticks_Y


	if xvisible == True:
		if ticks_X is not None:
			labelLon = ticks_Y
		else:
			labelLon = np.linspace(cord[2], cord[3], numTicks_Y)
			
		ax.set_xticks(labelLon, crs=proj)
		lon_formatter = LongitudeFormatter(
				zero_direction_label=zero_direction_label,
				dateline_direction_label=dateline_direction_label,
				number_format=number_format, degree_symbol=degree_symbol)
		ax.xaxis.set_major_formatter(lon_formatter)
		
	if yvisible == True:
		if ticks_Y is not None:
			labelLat = ticks_X
		else:
			labelLat = np.linspace(cord[0], cord[1], numTicks_X)
		ax.set_yticks(labelLat, crs=proj)
		lat_formatter = LatitudeFormatter(number_format=number_format,
				degree_symbol=degree_symbol)		
		ax.yaxis.set_major_formatter(lat_formatter)

	# set the label fontsize
	ax.tick_params( labelsize = fontsize)
	
	if y_label is not None:
		ax.set_ylabel(y_label, fontsize=fontsize)
	if x_label is not None:
		ax.set_xlabel(x_label, fontsize=fontsize)

		
	return ax

#-----------------------------------------------------------------------
def h_1_ax(fig, ax, x_off=0.0, y_off=-0.1, ratio=0.06,
           width_add=0.0, cbarWidth = None, cbarHeight=None,
           orientation = 'horizontal'):
	""" Add a colobar axes that is under one axes
	for *fig*.

	Parameters
	----------
	fig : plt.figure() instance
	ax : axes
	x_off : float (default is -0.1)
	   Offset from *ax* left
	y_off : float (default is -0.1)
	   Offset from *ax* bottom
	ratio : float (default is 0.06)
	   The ratio of height to width of the colorbar axes when
	   the *height* is not specified
	width_add : float (defaut is 0.0)
	   Default width colorbar axes equals that of *ax*.
	   *width_add* is used to change the width.
	height : float (default is None)
	   The height of the colorbar axes, it override *ratio* when
	   *height* is not None, or *height* is determined by *ratio*

	Returns
	 -------
	cax : axes
		A colobar axes that is under one axes for *fig*

	"""

	# axes poaistion
	pos = ax.get_position()

	width = pos.width  + width_add
	height = pos.height + width_add
	if orientation == 'horizontal':
		if cbarHeight is None:
			cbarHeight = width * ratio
		cax = fig.add_axes([pos.x0 , pos.y0 + y_off, width, cbarHeight])

	if orientation == 'vertical':
		if cbarWidth is None:
			cbarWidth = width * ratio
		cax = fig.add_axes([pos.x1 + x_off, pos.y0 , cbarWidth, height])	
# 	print('h_1_ax', pos.x0, pos.x1, pos.y0, pos.y1)
	return cax	
	
#-----------------------------------------------------------------------
def h_2_ax(fig, ax1, ax2, y_off=-0.1, ratio=0.06, height=None):
	""" Add a colobar axes that is centerd under two axes
	for *fig*.

	Parameters
	----------
	fig : plt.figure() instance
	ax1 : axes
	   left axes
	ax2 : axes
	   right axes
	y_off : float (default is -0.1)
	   Offset from ax1(ax2) bottom
	ratio : float (default is 0.06)
	   The ratio of height to width of the colorbar axes when
	   the *height* is not specified
	height : float (default is None)
	   The height of the colorbar axes, it override *ratio* when
	   *height* is not None, or *height* is determined by *ratio*

	Returns
	-------
	cax : axes
		A colobar axes that is centerd under two axes
		for *fig*

	"""

	# axes poaistion
	pos1 = ax1.get_position()
	pos2 = ax2.get_position()

	# cax
	x0 = pos1.x0 + pos1.width / 2.0
	y0 = pos1.y0 + y_off
	width = pos2.x0+pos2.width / 2.0 - (pos1.x0 + pos1.width / 2.0)
	if height is None:
		height = width * ratio
	cax = fig.add_axes([x0, y0, width, height])

	return cax

#-----------------------------------------------------------------------
def right_center_label(ax, label):
	""" Add *label* to the right center of *ax*.
	This function is usually used to add unit label
	on a colorbar

	Parameters
	----------
	ax : axes
	label : str

	"""

	ax.yaxis.set_label_position('right')
	ax.set_ylabel(label, rotation=0, ha='left', va='center')
		
#-----------------------------------------------------------------------
def plot_shapefile(ax, shapeData, **kwargs):

	c	= kwargs.get('c', 'orchid')
	ls  = kwargs.get('ls', 'dashdot')
	lw  = kwargs.get('lw', 'dashdot')
	aplha = kwargs.get('alpha', 0.5)
	
	for latList, lonList in zip(shapeData['latitude'], shapeData['longitude']):
		for i in range(len(latList)):
			lat = latList[i]
			lon = lonList[i]
			ax.plot(lon, lat, linewidth = lw, c = c, \
			        ls = ls, alpha =aplha)	
	return ax

#-----------------------------------------------------------------------
def explore_shapefile(shapeFile, **kwargs):

	import cartopy.io.shapereader as shpreader
	
	level = kwargs.get('cord', 'other')
	
	outlist = []
	
	if level == 'Detail':
		print('Detail')
	else:
		shp = shpreader.Reader(shapeFile)
		geometries = shp.geometries()
		countries = shp.records()
		for record, geometry in zip(countries, geometries):
			printSTR = ', '.join(record.attributes.keys())
		
			print(' - Explore_shapefile: Metadata of the shapefile, ', printSTR)
			for key in record.attributes.keys():
				outlist.append(key)
			break

		return outlist

#-----------------------------------------------------------------------
def read_shape(shapeFile, vars, **kwargs):

	import cartopy.io.shapereader as shpreader
	
	type = kwargs.get('type', 'county')
	
	shp = shpreader.Reader(shapeFile)
	geometries = shp.geometries()
	countries = shp.records()

	shapeData = {}
	
	# get the names of the attributes in the shapefile
	varList = explore_shapefile(shapeFile)
	
	# screen out the invalided variables...
	valid_var = []
	shapeData = {}
	for var in vars:
		if var in varList:
			valid_var.append(var)
			shapeData[var] = []
		else:
			print(' - Read_shape: ', var, 'is not in the shapefile, SKIP!')
	
	# create the a list to store the coordinates	
	shapeData['latitude'] = []
	shapeData['longitude'] = []
	j = 0
	for record, geometry in zip(countries, geometries):
		j = j + 1
		for var in valid_var:
			shapeData[var].append(record.attributes[var])
		latList = []
		lonList = []
		for i in range(len(geometry)):
			cord_poly = geometry[i]
			if type == 'county':
				# read the coordinates of polygen
				lons, lats = cord_poly.exterior.coords.xy
			elif type == 'street':
				# read the coordinates of linestring
				lons, lats = cord_poly.coords.xy
			else:
				print(' - Read_shape: type of coordinates does not exist...')
				lats = np.nan
				lons = np.nan		
			latList.append(lats)
			lonList.append(lons)
		shapeData['latitude'].append(latList)
		shapeData['longitude'].append(lonList)
	print(' - Read_shape: ', j, 'coordinates are found...')
		
	return shapeData

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
def plot_nasa_truecolor(cord, date_str, **kwargs):
	
	'''
	Function to plot the NASA world view true color image on a map. The default 
	requested turecolor image from NASA is 
	
	layer = VIIRS_SNPP_CorrectedReflectance_TrueColor
	other option is:
	
	1. MODIS_Terra_SurfaceReflectance_Band143
	2. MODIS_Terra_SurfaceReflectance_TrueColor
	3. MODIS_Aqua_SurfaceReflectance_Band143
	4. MODIS_Aqua_SurfaceReflectance_TrueColor
	
	Other false color images are also available, check them at 	
	http://esmc.uiowa.edu:3838/GIBS-leaflet-master/
	
	Parameters:
	cord     : list of coordinate of the image, [north, south, west, east]
	date_str : string, date of the requested image in a format of yyyy-mm-dd
	
	Returns:
	fig		: figure instance defined in this function
	ax		: axes instance defined in thes function
	plot_CRS: projection method defined in thes function
	
	'''
	from owslib.wmts import WebMapTileService
	import numpy as np
	import cartopy.crs as ccrs
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
	import cartopy.feature as cfeature	
	
	ysize = kwargs.get('figsize', 8)
	title		= kwargs.get('title', '')
	
	numTicks_Y  = kwargs.get('numTicks_Y', 5)
	numTicks_X  = kwargs.get('numTicks_X', 5)
	
	linewidth  	= kwargs.get('linewidth', 0.75)
	fontsize  	= kwargs.get('fontsize', 12)
	lineColor	= kwargs.get('lineColor', 'lightgray')
	dateline_direction_label = kwargs.get('dateline_direction_label', False) 
	degree_symbol = kwargs.get('degree_symbol', u'\u00B0') 
	number_format = kwargs.get('number_format', '4.0f') 
	zero_direction_label = kwargs.get('zero_direction_label', False)
	cbar_pos = kwargs.get('cbar_pos', 'vertical') 
	yvisible = kwargs.get('yvisible', True)
	xvisible = kwargs.get('xvisible', True)
	
	# Layers for nasa gibs...
	layer = kwargs.get('layer', 'VIIRS_SNPP_CorrectedReflectance_TrueColor')	
	
	
	URL = 'https://gibs.earthdata.nasa.gov/wmts/epsg4326/best/wmts.cgi'
	wmts = WebMapTileService(URL)

	
	# Plot setup
	plot_CRS = ccrs.PlateCarree(central_longitude=0)
	geodetic_CRS = ccrs.Geodetic()
	x0, y0 = plot_CRS.transform_point(cord[2], cord[1], geodetic_CRS)
	x1, y1 = plot_CRS.transform_point(cord[3], cord[0], geodetic_CRS)


	labelLat = np.linspace(cord[1], cord[0], numTicks_Y)
	labelLon = np.linspace(cord[2], cord[3], numTicks_X)

	# calculate the x size base on the ratio of the domain size
	xsize = 2 * ysize * (x1 - x0) / (y1 - y0)
	axsize = [0.15,0.15,0.7,0.7]

	fig = plt.figure(figsize=(xsize, ysize),  dpi=200)
	ax = fig.add_axes(axsize, projection=plot_CRS)
	ax.set_xlim((x0, x1))
	ax.set_ylim((y0, y1))
	ax.add_wmts(wmts, layer, wmts_kwargs={'time': date_str})

	# set x ticks
	if xvisible == True:
		ax.set_xticks(labelLon, crs=plot_CRS)
		lon_formatter = LongitudeFormatter(
				zero_direction_label=zero_direction_label,
				dateline_direction_label=dateline_direction_label,
				number_format=number_format, degree_symbol=degree_symbol)
		ax.xaxis.set_major_formatter(lon_formatter)

	# set y ticks
	if yvisible == True:
		ax.set_yticks(labelLat, crs=plot_CRS)
		lat_formatter = LatitudeFormatter(number_format=number_format,
				degree_symbol=degree_symbol)		
		ax.yaxis.set_major_formatter(lat_formatter)

	# add the states boundary here
	ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=linewidth, edgecolor = lineColor)

	# add the lake boundary here
	ax.add_feature(cfeature.LAKES.with_scale('10m'), linewidth=linewidth, \
				   edgecolor = lineColor, alpha=0.5)
	# add the river line here
	ax.add_feature(cfeature.RIVERS.with_scale('10m'), linewidth=linewidth, \
				   edgecolor = lineColor)
			   
	# add the coastline, including major islands.
	ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=linewidth, \
				   edgecolor = lineColor)
	# add the land polygons, including major islands.
	ax.add_feature(cfeature.LAND.with_scale('10m'), linewidth=linewidth, \
				   edgecolor = lineColor)
	# add the country boundaries.
	ax.add_feature(cfeature.BORDERS.with_scale('10m'), linewidth=linewidth, \
				   edgecolor = lineColor,linestyle=':')
				   
	ax.tick_params( labelsize = fontsize)			   
	# set the title of the figure
	ax.set_title(title, fontsize = fontsize)
	
	
	return fig, ax, plot_CRS


#-----------------------------------------------------------------------
def plotSat_2d_backgroud( ax, proj, cord, **kwargs):

	'''
	Function to create the basic basemap projection figure...
    
	Parameters
    ----------    
	ax		: axis instance
	
	proj	: basemap instance      
    cord	: coordinates of the map, [north, south, west, east]

	Return
    ----------  
	ax		: axis instance

	optional Parameters:
	----------  
	cord		: list	
	scaleFactor	: float
	vRange		: list
	cbar_label  : str
	title		: str
	numTicks_Y  : int
	numTicks_X  : int
 	lineWidth  	: float
 	fontsize  	: float
 	lineColor	: str
 	cmap		: str
 	badColor	: str
	
	'''

	import numpy as np
	import matplotlib.pyplot as plt
	import cartopy.crs as ccrs
	import cartopy.feature as cfeature
	from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
	
	
	# default parameters are set here	
	title		= kwargs.get('title', '')
	numTicks_Y  = kwargs.get('numTicks_Y', 5)
	numTicks_X  = kwargs.get('numTicks_X', 5)
	linewidth  	= kwargs.get('linewidth', 0.6)
	fontsize  	= kwargs.get('fontsize', 14)
	lineColor	= kwargs.get('lineColor', 'dimgray')
	numColr		= kwargs.get('numColr', 255)	
	dateline_direction_label = kwargs.get('dateline_direction_label', False) 
	degree_symbol = kwargs.get('degree_symbol', u'\u00B0') 
	number_format = kwargs.get('number_format', '4.2f') 
	zero_direction_label = kwargs.get('zero_direction_label', False)
	cbar_pos = kwargs.get('cbar_pos', 'vertical') 
	yvisible = kwargs.get('yvisible', True)
	xvisible = kwargs.get('xvisible', True)

	# create the label of the x and y axis
	labelLat = np.linspace(cord[0], cord[1], numTicks_Y)
	labelLon = np.linspace(cord[2], cord[3], numTicks_X)
	
	# set the paralle line of the latitude and longitude, I personal prefer the line 
	# width is half width of other line...	
	
	if xvisible == True:

		ax.set_xticks(labelLon, crs=proj)
		lon_formatter = LongitudeFormatter(
				zero_direction_label=zero_direction_label,
				dateline_direction_label=dateline_direction_label,
				number_format=number_format, degree_symbol=degree_symbol)
		ax.xaxis.set_major_formatter(lon_formatter)
		
	if yvisible == True:
		ax.set_yticks(labelLat, crs=proj)
		lat_formatter = LatitudeFormatter(number_format=number_format,
				degree_symbol=degree_symbol)		
		ax.yaxis.set_major_formatter(lat_formatter)

	ax.gridlines(linewidth=linewidth*0.75, crs=proj, xlocs=labelLon, ylocs=labelLat, 
	             edgecolor = lineColor,linestyle=':')

	
	ax.tick_params( labelsize = fontsize)

	# add the states boundary here
	ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=linewidth, edgecolor = lineColor)

	# set the title of the figure
	ax.set_title(title, fontsize = fontsize)

	return ax

#-----------------------------------------------------------------------
def calPixelEdge(data):
	"""
	using center latitude and longitude to calculate edge
	"""
	import numpy as np
	dim = data.shape
	num_l = dim[0]
	num_p = dim[1]

	data_e = np.zeros((num_l+2,num_p+2))

	data_e[1:-1,1:-1] = data


	# extend the left column...
	data_e[1:-1,0]   = 2*data[:,0] - data[:,1]
	# extend the right column
	data_e[1:-1,-1]   = 2*data[:,-1] - data[:,-2]
	
	data_e[0,:]       = 2*data_e[1,:] - data_e[2,:]
	data_e[-1,:]   = 2*data_e[-2,:] - data_e[-3,:]
	data_e = (data_e[0:-1, 0:-1] + data_e[1:, 1:] + data_e[1:, 0:-1] + data_e[0:-1, 1:])/4.
# 	print(data_e)
# 	
# 	data_e = np.zeros((num_l+1,num_p+1))
# 	
# 	data_e_tmp = (data[0:num_l-1,:] + data[1:num_l,:]) * 0.5
# # 	print('1', data_e_tmp)
# 	data_e_tmp = (data_e_tmp[:,0:num_p-1] + data_e_tmp[:,1:num_p]) * 0.5
# # 	print('2',data_e_tmp)
# 	data_e[1:num_l,1:num_p] = np.array(data_e_tmp)
# # 	print('3',data_e)
# 	data_e[0,:]  = 2.0 * data_e[1,:]  - data_e[2,:]
# 	data_e[-1,:] = 2.0 * data_e[-2,:] - data_e[-3,:]
# 	data_e[:,0]  = 2.0 * data_e[:,1]  - data_e[:,2]
# 	data_e[:,-1] = 2.0 * data_e[:,-2] - data_e[:,-3]

	return data_e

#-----------------------------------------------------------------------
def plot_patch(ax, GridCell, **kwargs):
	
	'''
	Function to plot a bunch of grids on a map
	-----------
	Parameters:
	
	ax		 :	axes instance
	
	GridCell : coordinates of the grids
	
	-----------
	Returns:
	 
	ax		 :	axes instance
	
	Key optional parameters:
	
	faceColorOn	: if Ture, front color of the pixel will be determined by
				  colorValue, otherwise it will fill the all the pixels 
				  the same color
	cmap		: c
	'''
	
	import matplotlib.cm as cm
	import matplotlib.pyplot as plt
	from matplotlib.path import Path
	import matplotlib.colors as colors
	import matplotlib.patches as patches
	from matplotlib.collections import PatchCollection
	

	# default parameters are set here
	vRange		= kwargs.get('vRange', [0,1])
	edgecolors	= kwargs.get('edgecolors', 'black')
	linewidth	= kwargs.get('linewidth', '2')
	faceColorOn = kwargs.get('faceColorOn', False)
	colorValue 	= kwargs.get('colorValue', None)
	cmap		= kwargs.get('cmap', None)
	cLognorm	= kwargs.get('cLognorm', False)
	numColr		= kwargs.get('numColr', 255)
	x_idx		= kwargs.get('x_idx', None)
	y_idx		= kwargs.get('y_idx', None)
	alpha		= kwargs.get('alpha', 1.0)
	linestyle	= kwargs.get('linestyle', '-')
	fill		= kwargs.get('fill', True)
	
	# - get the coordinates of the pixel edges
	if x_idx == None:
		ullon= GridCell['Upper_Left_Lon']
		urlon= GridCell['Upper_Right_Lon']
		lllon= GridCell['Lower_Left_Lon']
		lrlon= GridCell['Lower_Right_Lon']

		ullat= GridCell['Upper_Left_Lat']
		urlat= GridCell['Upper_Right_Lat']
		lllat= GridCell['Lower_Left_Lat']
		lrlat= GridCell['Lower_Right_Lat']
	else:
		ullon= GridCell['Upper_Left_Lon'][x_idx, y_idx].flatten()
		urlon= GridCell['Upper_Right_Lon'][x_idx, y_idx].flatten()
		lllon= GridCell['Lower_Left_Lon'][x_idx, y_idx].flatten()
		lrlon= GridCell['Lower_Right_Lon'][x_idx, y_idx].flatten()

		ullat= GridCell['Upper_Left_Lat'][x_idx, y_idx].flatten()
		urlat= GridCell['Upper_Right_Lat'][x_idx, y_idx].flatten()
		lllat= GridCell['Lower_Left_Lat'][x_idx, y_idx].flatten()
		lrlat= GridCell['Lower_Right_Lat'][x_idx, y_idx].flatten()

	patchess = []
	Marea = []
	
	# prepare the color setting for each pixel
	if faceColorOn == False:
		# set all the pixel the same color
		colorValue = []
		for i in range(len(ullon)):
			colorValue.append(cmap)
	else:
	 	# otherwise set color based on pixel value
		cmap = plt.get_cmap(cmap, numColr)
		vmin = vRange[0]
		vmax = vRange[1]
		if cLognorm:
			# log scale color...
			cNorm  = colors.LogNorm(vmin=vmin, vmax=vmax)
			scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
		
		else:
			# otherwise, normally distribute the color
			cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
			scalarMap = cm.ScalarMappable(norm=cNorm, cmap=cmap)
			
	# create the pixel instance
	for i,j,k,m,n,o,p,s,clr in zip(ullon,urlon,lllon,lrlon,\
	                               ullat,urlat,lllat,lrlat,\
	                               colorValue):
		verts = [
				(k, p),		# left, bottom
				(i, n),		# left, top
				(j, o), 	# right, top
				(m, s), 	# right, bottom
				(0., 0.), 	# ignored
				]
		codes = [Path.MOVETO, Path.LINETO, Path.LINETO, \
				 Path.LINETO, Path.CLOSEPOLY,]
		path = Path(verts, codes)
		
		# create the pixel
		if faceColorOn:

			facecolor = scalarMap.to_rgba(clr)
			rectangle = patches.PathPatch(path, facecolor=facecolor, \
										  linewidth=linewidth,  \
										  edgecolor=edgecolors, \
										  alpha = alpha,\
										  linestyle=linestyle)
		else:

			facecolor = clr
			scalarMap = clr

			if facecolor is None:
				rectangle = patches.PathPatch(path, facecolor=None, \
											  linewidth=linewidth, \
											  edgecolor = edgecolors, \
											  fill = False, \
											  alpha = alpha,\
											  linestyle=linestyle)
			else:
				rectangle = patches.PathPatch(path, facecolor=facecolor, \
											  linewidth=linewidth, \
											  edgecolor = edgecolors, \
											  fill=fill, \
											  alpha = alpha,\
											  linestyle=linestyle)

		# add the patch into the axes
		ax.add_patch(rectangle)
		patchess.append(rectangle)

	return ax, scalarMap

#-----------------------------------------------------------------------
def set_colorbar(fig, ax, cbar_pos, **kwargs):

	'''
	Function to add the color bar instance on a figure...
	
	------------
	Parameters:
	fig	 	 : figure instance 
	ax	 	 : axes instance  
	cbar_pos : orientation of color bar, vertical or horizontal
	
	Other optional parameters are explained in the code
	
	------------
	Returns:
	ax	 	 : axes instance  
	
	'''
	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
		
	# related to the size and position of the colorbar...
	x_off = kwargs.get('x_off', 0.015)
	y_off = kwargs.get('y_off', 0)
	cbarRatio = kwargs.get('cbarRatio', 0.05)	
	# related to the position of the color bard label
	labelpad 	= kwargs.get('labelpad', 1.5)
	# related to the extend of the color bar, max/min/neither/both
	extend		= kwargs.get('extend', 'neither')
	# related to the fontsize of the color bar
	fontsize = kwargs.get('fontsize', 12)
	
	# label of the color bar
	cbar_label = kwargs.get('cbar_label', '')
	
	# number of the color used to denote the value
	numClr = kwargs.get('numClr', 255)
	# number of the ticks of the color bar
	numCTick = kwargs.get('numCTick', 5)
	# range of value of the color bar
	vRange = kwargs.get('vRange', [0, 1])
	
	# logical switch, if True, the color bar will denote as log-scale,
	# this term will only modify the color bar ticks, it won't change 
	# actual the value of the actual dataset... 
	cLognorm = kwargs.get('cLognorm', False)
	
	# format of the the color bar ticks
	cbarFomat = kwargs.get('cbarFomat', r'%4.2f')
	
	# list-like, if provided, the ticks of the color bar will 
	# determined by this term
	
	cbarTick       = kwargs.get('cbarTick', None)
	
	cbarTickLabels = kwargs.get('cbarTickLabels', None)
	
	# get the color range (maximum, minimum) of the data
	vmin = vRange[0]
	vmax = vRange[1]

	# get the position of the color bar	
	cax_v = h_1_ax(fig, ax, x_off = x_off, y_off= y_off, orientation = cbar_pos, ratio = cbarRatio)
	# create the color bar instance
	cbar = plt.colorbar(ax.collections[0], cax=cax_v, orientation = cbar_pos, extend=extend)
	
	# set the fontsize of the ticks
	cbar.ax.tick_params(labelsize = fontsize)
	# set the label
	cbar.set_label(cbar_label, fontsize = fontsize, labelpad=labelpad)
	
	# two options of the ticks of the color bar
	# I.  discrete color
	# II. continuous color
	if numClr < 255:
		# a discrete color bar is selected...
		if cbarTick == None:
			deltaTick = np.abs(vmax - vmin)/numCTick/2
			cbarTick  = np.linspace(vmin, vmax, numCTick+1) + deltaTick
			cbarTick  = cbarTick[0:numCTick]


		if cbarTickLabels == None:
			# if the color bar ticks is not provided...then generate the 
			# ticks based on the range of the value, and number of the ticks
			cbarTickLabels = [cbarFomat%tick for tick in np.linspace(vmin, vmax, numCTick)]
		
		# set the tick in the middle of the color block
		cbar.set_ticks(cbarTick)
		cbar.set_ticklabels(cbarTickLabels)
	else:
		# a continuous color bar is selected...
		if cLognorm:
			# denote the color in a log-10 scale...
			vmin = np.log10(vRange[0])
			vmax = np.log10(vRange[1])
			if cbarTick == None:
				cbarTick = np.logspace(vmin, vmax, numCTick)
			if cbarTickLabels == None:
				cbarTickLabels = [cbarFomat%tick for tick in cbarTick]
		else:
			# denote the value in normal scale...
			if cbarTick is None:
				cbarTick = np.linspace(vmin, vmax, numCTick)
			# denote the value in normal scale...	
			if cbarTickLabels == None:
				# if color bar ticks is not provided...then generate the ticks based
				# on the range of the value, and number of the ticks...				
				cbarTickLabels = [cbarFomat%tick for tick in cbarTick]
	
		cbar.set_ticks(cbarTick)
		cbar.set_ticklabels(cbarTickLabels)
			
	return ax


#-----------------------------------------------------------------------

def patchOnMap(fig, ax, proj, GridCells, **kwargs):

	colorValue		= kwargs.get('colorValue', None)
	xvisible		= kwargs.get('xvisible', False)
	yvisible		= kwargs.get('yvisible', False)
	numCTick		= kwargs.get('numCTick', 4)
	cbar_label		= kwargs.get('cbar_label', '')
	colorExtend		= kwargs.get('colorExtend', 'both')
	cbarRatio		= kwargs.get('cbarRatio', 0.05)
	cbar_pos		= kwargs.get('cbar_pos', 'vertical')
	cbarOn			= kwargs.get('cbarOn', False)
	cmap			= kwargs.get('cmap', 'rainbow')
	cLognorm		= kwargs.get('cLognorm', False)
	vRange			= kwargs.get('vRange', [0,1])
	faceColorOn		= kwargs.get('faceColorOn', False)
	numTicks_X		= kwargs.get('numTicks_X', 5)
	numTicks_Y		= kwargs.get('numTicks_Y', 5)
	number_format 	= kwargs.get('number_format', '')	
	fontsize		= kwargs.get('fontsize', '4.2f')
	level			= kwargs.get('level', '5')
	title			= kwargs.get('title', '')
	cbarFomat		= kwargs.get('cbarFomat', r'%4.3f')
	numColr			= kwargs.get('numColr', 255)
	gisOn			= kwargs.get('gisOn', False)
	cord			= kwargs.get('cord', [90, -90, -180, 180])

	if gisOn:
		ax = plot_2d_gis(ax, proj, level, cord = cord, **kwargs )
					 
	ax, scalarMap = plot_patch(ax, GridCells, **kwargs)
	if cbarOn:
		ax = set_colorbar(fig, ax, scalarMap, cbar_pos, **kwargs)
	
	return ax

#-----------------------------------------------------------------------
def customize_colormap(rgb_seq, num_color, **kwargs):
	"""
	customize a colormap
	
	------------
	Parameters:
	rgb_seq : list of list, RGB color, 
	          rgb_seq = [[  0,  0,200],	   
		                 [255,255,255], # white
                         [217,  3, 24]] 
    num_color: number of color that will be used in the colormap
	
	------------
	Returns:
	data_e: array, edge of the input latitude or longitude
	
	"""

	import numpy as np
	from matplotlib.colors import LinearSegmentedColormap
	# if Ture, will create an uneven color
	uneven    = kwargs.get('uneven', False) 
	# number of colors that will be use below the center color...
	under_Num = kwargs.get('under_Num', 5)

	
	if uneven == False:
		cm_data = []

		for colors in rgb_seq:
			uni_color = []
			for color in colors:
				uni_color.append(color/255.)
			cm_data.append(uni_color)	
		cm = LinearSegmentedColormap.from_list('my_color', cm_data, num_color)
		
	if uneven == True:
		
		over_Num = num_color - under_Num 
		cm_data = []
		
		
		for colors in rgb_seq:
			uni_color = []
			for color in colors:
				uni_color.append(color/255.)
			cm_data.append(uni_color)	
		cm_1 = LinearSegmentedColormap.from_list('my_color', cm_data, 255)		
		
		
		cm_data = []
		
		color_pos = np.linspace(0,0.5, under_Num + 1 )
		
		for i in range(under_Num):
			cm_data.append( cm_1(color_pos[i]) )
			
		color_pos = np.linspace( 0.5, 1.0, over_Num + 1)	
		cm_data.append( cm_1(0.5) )	
		for i in range(over_Num):		
			cm_data.append( cm_1(color_pos[i+1]) )
		cm = LinearSegmentedColormap.from_list('my_color', cm_data, under_Num + over_Num + 1)	
		
	return cm

#-----------------------------------------------------------------------
def bytescale(data, cmin=None, cmax=None, high=255, low=0):
	"""
	Byte scales an array (image).

	Byte scaling means converting the input image to uint8 dtype and 
	scaling the range to ``(low, high)`` (default 0-255).
	If the input image already has dtype uint8, no scaling is done.

	Parameters
	----------
	data : ndarray
		PIL image data array.
	cmin : scalar, optional
		Bias scaling of small values. Default is ``data.min()``.
	cmax : scalar, optional
		Bias scaling of large values. Default is ``data.max()``.
	high : scalar, optional
		Scale max value to `high`.  Default is 255.
	low : scalar, optional
		Scale min value to `low`.  Default is 0.

	Returns
	-------
	img_array : uint8 ndarray
		The byte-scaled array.


	"""
	import numpy as np
	if data.dtype == np.uint8:
		return data

	if high < low:
		raise ValueError("`high` should be larger than `low`.")

	if cmin is None:
		cmin = np.nanmin(data)
	if cmax is None:
		cmax = np.nanmax(data)

	cscale = cmax - cmin
	if cscale < 0:
		raise ValueError("`cmax` should be larger than `cmin`.")
	elif cscale == 0:
		cscale = 1

	scale = float(high - low) / cscale
	bytedata = (data * 1.0 - cmin) * scale + 0.4999
	# set the NaN value to 0 to avoid warning
	bytedata[data!=data] = 0
	bytedata[bytedata > high] = high
	bytedata[bytedata < 0] = 0
	bytedata = np.cast[np.uint8](bytedata) + np.cast[np.uint8](low)
	bytedata = np.ma.masked_array(bytedata, data!=data)
# 	bytedata[data!=data] = 0
	return bytedata

#-----------------------------------------------------------------------
def bytescale_multi_figure(data_list, channels):
	import numpy as np
	cmin = 1000
	cmax = 0
	print('\n - Calling bytescale_multi_figure')
	channel_max = {}
	channel_min = {}
	
	for channel in channels:
		channel_max[channel] = cmax
		channel_min[channel] = cmin
	
	
	for data in data_list:
		for channel in channels:
			cmin = np.nanmin(data[channel])
			cmax = np.nanmax(data[channel])
			
			if channel_max[channel] < cmax:
				channel_max[channel] = cmax
			if channel_min[channel] > cmin:
				channel_min[channel] = cmin
	

	for data in data_list:
		for channel in channels:
			cmax = channel_max[channel]
			cmin = channel_min[channel]
			print('\n - Channel', channel, 'max:%6.2f, mmin:%6.2f'%(cmax, cmin) )
			data[channel] = bytescale(data[channel], cmin=cmin, cmax=cmax)
	
	return data_list	
#-----------------------------------------------------------------------
def scale_image(image, x=None, y=None):
	""" 	
	-----------
	Parameters:
	image: 2D array, image
	x	 : the originate brightness degree of the image
	y	 : the target brightness degree of the image 
	
	-----------
	Returns:
	image: 2D array, image
			
	color enhancement
	https://moonbooks.org/Codes/Plot-MODIS-granule-RGB-image-
	orthographic-projection--color-enhancement-using-python-and-basemap/
	
	"""
	import numpy as np 
	if x is None:
		x = np.array([0,  30,  60, 120, 190, 255], dtype=np.uint8)

	if y is None:
		y = np.array([0, 110, 160, 210, 240, 255], dtype=np.uint8)

	# bytescale
	image_bs = bytescale(image)

	scaled = np.zeros_like(image_bs, dtype=np.uint8)
	for i in range(len(x)-1):
		x1 = x[i]
		x2 = x[i+1]
		y1 = y[i]
		y2 = y[i+1]
		m = (y2 - y1) / float(x2 - x1)
		b = y2 - (m *x2)
		mask = ((image_bs >= x1) & (image_bs < x2))
		scaled = scaled + mask * np.asarray(m * image_bs + b, dtype=np.uint8)

	mask = image_bs >= x2
	scaled = scaled + (mask * 255)

	scaled = scaled.astype(float)

	scaled = scaled / 255.0

	return scaled
#-----------------------------------------------------------------------	
def scale_image_multi_c(image, x=None, y=None):
	""" 
	color enhancement of mutliple channnels.
	"""
	import numpy as np
	print('\n - Calling scale_image_multi_c')
	scaled = np.zeros_like(image)

	for i in range(image.shape[-1]):
		scaled[:,:,i] = scale_image(image[:,:,i], x=x, y=y)

	return scaled
	
	
#-----------------------------------------------------------------------
def convolution2d(image, **kwargs):

	filter_type		= kwargs.get('filter_type', None)
	
	import numpy as np
	'''
	Function to conduct a 2D convolution using the kernel provided
	
	
	----------
	parameters:
	image   : 2D array of image
	kernel  : 2D template operator
	
	----------
	return 
	re_img  : convolution results

	'''
	m, n = 3,3

# 	m, n = kernel.shape
	y, x = image.shape
	
	l = (m-1)//2
	new_y = y + 2*l 
	new_x = x + 2*l 
	
	new_image = np.zeros((new_y,new_x))
	new_image[l :-l , l :-l ] = image
	
	re_img = np.zeros((y, x))
	
	for i in range(l,new_y-l):
		for j in range(l,new_x-l):
			if filter_type == None:
				re_img[i-l][j-l] = np.nanmean(new_image[i-l:i+l+1, j-l:j+l+1] )
			if filter_type == 'min':
				re_img[i-l][j-l] = np.nanmin(new_image[i-l:i+l+1, j-l:j+l+1] )
			if filter_type == 'med':
				re_img[i-l][j-l] = np.nanmedian(new_image[i-l:i+l+1, j-l:j+l+1] )
	return re_img