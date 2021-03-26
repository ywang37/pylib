# -*- coding: utf-8 -*-
"""
Modified Taylor diagram (Taylor, 2001) test implementation.
see 

Created on Fri Dec 05 15:39:18 2014

@author: laurence
inspired by:
"Yannick Copin's <yannick.copin@laposte.net> Taylor Diagram implimentation
  https://gist.github.com/ycopin/3342888"
  
Sean Elvidge's European Space Weather Week 11 (2014) talk 
  http://www.stce.be/esww11/contributions/public/splinters/metrics/SeanElvidge/

Elvidge, S., M. J. Angling, and B. Nava (2014), 
  On the use of modified Taylor diagrams to compare ionospheric assimilation models, 
  Radio Sci., 49, 737–745, doi:10.1002/2014RS005435.
  
Taylor, K.E. (2001):  
  Summarizing multiple aspects of model performance in a single diagram. 
  J. Geophys. Res., 106, 7183-7192 
  (also see PCMDI Report 55, http://www-pcmdi.llnl.gov/publications/ab55.html)

see also
http://www-pcmdi.llnl.gov/about/staff/Taylor/CV/Taylor_diagram_primer.htm  

Modified BSD License

Copyright (c) 2014 Laurence Billingham.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

  a. Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
  b. Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
  c. Neither the name of the Scikit-learn Developers  nor the names of
     its contributors may be used to endorse or promote products
     derived from this software without specific prior written
     permission. 


THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
POSSIBILITY OF SUCH DAMAGE.
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist.angle_helper as angle_helper
import mpl_toolkits.axisartist.floating_axes as fa
import mpl_toolkits.axisartist.grid_finder as gf
from matplotlib import colors as mcolors
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import SubplotHost
from mpl_toolkits.axisartist import GridHelperCurveLinear

class TaylorDiagramPoint(object):
  """
  A single point on a Modified Taylor Diagram.
  How well do the values predicted match the values expected
      
      * do the means match
      * do the standard deviations match
      * are they correlated
      * what is the normalized error standard deviation
      * what is the bias?
  
  Notation:
  
      * s_ := sample standard deviation 
      * m_ := sample mean 
      * nesd := normalized error standard deviation;
                > nesd**2 = s_predicted**2 + s_expected**2 -
                            2 * s_predicted * s_expected * corcoeff
  """
  def __init__(self, expected, predicted, pred_name, point_id, marker_id,
          normalized):

    self.normalized = normalized
    self.pred = predicted
    self.expd = expected
    self.s_pred = np.std(self.pred)
    self.s_expd = np.std(self.expd)
    if self.normalized:
        self.s_normd = self.s_pred / self.s_expd
    else:
        self.s_normd = self.s_pred
    self.bias = (np.mean(self.pred) - np.mean(self.expd)) / np.mean(self.expd) * 100
    self.corrcoef = np.corrcoef(self.pred, self.expd)[0, 1]
    self.corrcoef = min([self.corrcoef, 1.0])
#    self.nesd = np.sqrt(self.s_pred**2 + self.s_expd**2 - 
#                   2 * self.s_pred * self.s_expd * self.corrcoef)
    self.name = pred_name
    self.point_id = point_id
    self.marker_id = marker_id

class ModTaylorDiagram(object):
  """
    Given predictions and expected numerical data 
    plot the standard deviation of the differences and correlation between
    expected and predicted in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).
  """
  def __init__(self, fig=None, label='expected', max_normed_std=2.5, std_ratios=None, \
          bias_vmin=None, bias_vmax=None, normalized=True, \
          cmap=plt.cm.jet, s=80, title_expected=r'expected'):
    """
    Set up Taylor diagram axes. 
    """   
    
    self.normalized     = normalized
    self.title_polar    = r'Correlation'
    if self.normalized:
        self.title_xy   = r'Normalized Standard Deviation'
    else:
        self.title_xy   = r'Standard Deviation'
    self.max_normed_std = max_normed_std
    self.s_min          = 0
    self.std_ratios     = std_ratios
    self.bias_vmin      = bias_vmin
    self.bias_vmax      = bias_vmax
    self.cmap           = cmap
    self.s              = s # marker size
    self.title_expected = title_expected

    # Correlation labels
    corln_r = np.concatenate((np.arange(10)/10.,[0.95,0.99]))
    corln_ang = np.arccos(corln_r)      # Conversion to polar angles
    grid_loc1 = gf.FixedLocator(corln_ang)    # Positions
    tick_fmttr1 = gf.DictFormatter(dict(zip(corln_ang, map(str, corln_r))))
    
    # Normalized standard deviation axis
    tr = PolarAxes.PolarTransform()
    grid_helper = fa.GridHelperCurveLinear(tr,
                                  extremes=(0, np.pi/2, # 1st quadrant
                                            self.s_min, self.max_normed_std),
                                  grid_locator1=grid_loc1,
                                  tick_formatter1=tick_fmttr1)
    self.fig = fig                                   
    if self.fig is None:
      self.fig = plt.figure()
    
    # setup axes
    ax = fa.FloatingSubplot(self.fig, 111, grid_helper=grid_helper)
    # make the axis (polar ax child used for plotting)    
    self.ax = self.fig.add_subplot(ax)
    # hide base-axis labels etc   
    self.ax.axis['bottom'].set_visible(False)  
    self._setup_axes()
                   
    # attach the polar axes and plot the correlation lines
    self.polar_ax = self.ax.get_aux_axes(tr)
    self.polar_ax.plot([np.pi/2.135, np.pi/2.135], [self.s_min, self.max_normed_std], color='grey')
    self.polar_ax.plot([np.pi/2.484, np.pi/2.484], [self.s_min, self.max_normed_std], color='grey')
    self.polar_ax.plot([np.pi/3,     np.pi/3    ], [self.s_min, self.max_normed_std], color='grey')    
    self.polar_ax.plot([np.pi/3.95,  np.pi/3.95 ], [self.s_min, self.max_normed_std], color='grey')    
    self.polar_ax.plot([np.pi/6.95,  np.pi/6.95 ], [self.s_min, self.max_normed_std], color='grey')

    # Add norm stddev ratio contours
    if self.normalized:
        self._plot_req_cont()
    
    self.points = []

  def add_prediction(self, expected, predicted, predictor_name='', 
                     plot_pt_id='', plot_marker_id='o'):
    """
    Add a prediction/model to the diagram
    """
    this_point = TaylorDiagramPoint(expected, predicted, 
                                    predictor_name, plot_pt_id, plot_marker_id,
                                    self.normalized)   
    self.points.append(this_point)
    
  def plot(self):
    """
    Place all the loaded points onto the figure 
    """
    rs         = []
    std        = []
    thetas     = []
    biases     = []
    names      = []
    point_tags = []
    marker     = []
    for point in self.points:
      rs.append(point.s_expd)
      std.append(point.s_normd)
      thetas.append(np.arccos(point.corrcoef))
      biases.append(point.bias)
      names.append(point.name)
      point_tags.append(point.point_id)
      marker.append(point.marker_id)

    #Plot the angle values and normalized standard deviations
    # Add the point labels and offset them so they’re readable
    if self.bias_vmin is None:
        self.bias_vmin = min(biases)
    if self.bias_vmax is None:
        self.bias_vmax = max(biases)
    norm = mcolors.Normalize(vmin=self.bias_vmin,vmax=self.bias_vmax)
    for i in range(len(marker)):
        self.sc = self.polar_ax.scatter([thetas[i]], [std[i]], c=[biases[i]], \
                s=self.s, cmap=self.cmap, norm=norm, marker=marker[i], \
                edgecolors='k')
    for i, tag in enumerate(point_tags):                         
      self.polar_ax.text(thetas[i]+0.01, std[i]+0.01, tag,
                         horizontalalignment='center',
                         verticalalignment='bottom', fontsize=16)

              
  def show_key(self):
    """
    add annotation key for model IDs and normalization factors
    annotation key currently includes the point id (letter) and site name
    """
    textstr = ''
    for i, p in enumerate(self.points):
      if i > 0:
        textstr += '\n'
        
      textstr += r'{0}$\rightarrow${1}'.format(p.point_id, p.name)                        
    props = dict(boxstyle='round', facecolor='white', alpha=0.75)
    # place a text box in upper left in axes coords
    self.ax.text(1, 1, textstr, transform=self.ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

  def bias_colorbar(self, **kwargs):
      """
      add bias colorbar
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      cbar = plt.colorbar(self.sc, **kwargs)

      cbar.set_label('Relative bias [%]', fontsize=12)

      return cbar

  def stringID_legend(self, stringID_labels, **kwargs):
      """
      add stringID legend
      stringID_labels is tuple, and every element in the tuple is also a tuple which
      consists marker and label.
      for example, stringID_labels: ( (r'$1$', 'OMI'), (r'$2$', 'OMPS') )
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      handles = []
      for stID_la in stringID_labels:
          stID = stID_la[0]
          la   = stID_la[1]
          handles.append( plt.scatter([], [], marker=stID, label=la, s=self.s, c='k') )

      legend = plt.legend(handles=handles, **kwargs)
      self.ax.add_artist(legend)

  def marker_legend(self, marker_labels, **kwargs):
      """
      add marker legend:
      marker_labels is a tuple, and every element in the tuple is also a tuple which 
      consists marker and label.
      for example, marker_labels: ( ('o', 'prior'), ('s', 'posterior') )
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      handles = []
      for ma_la in marker_labels:
          ma = ma_la[0]
          la = ma_la[1]
          handles.append( plt.scatter([], [], marker=ma, label=la, edgecolors='k', s=self.s, c='w') )

      legend = plt.legend(handles=handles, **kwargs)
      self.ax.add_artist(legend)
        
  def show_norm_factor(self):
    """
    add annotation about the normalization factor
    """        
    n_fact = self.points[0]
    out_str = r'Norm Factor {:.2f}'.format(n_fact)   
    x = 0.95 * self.max_normed_std   
    y = 0.95 * self.max_normed_std
    self.ax.text(x, y, 
                 out_str,
                 horizontalalignment='right',
                 verticalalignment='top', 
                 bbox={'edgecolor': 'black', 'facecolor':'None'})      
                 
  def add_contours(self, levels=5, **kwargs):
    """
    Add constant centered RMS difference contours around reference point
    """

    rs,ts = np.meshgrid(np.linspace(self.s_min,self.max_normed_std), \
            np.linspace(0,np.pi/2))

    # Compute centered RMS difference
    if self.normalized:
        refstd = 1.0
    else:
        # not mormailized by expected standard deviation (s_expd)
        # (Yi Wang, 11/02/2019)
        refstd = self.points[0].s_expd
        rms = np.sqrt((refstd)**2 + rs**2 - 2*(refstd)*rs*np.cos(ts))
    rms = np.sqrt((refstd)**2 + rs**2 - 2*(refstd)*rs*np.cos(ts))
    
    contours = self.polar_ax.contour(ts, rs, rms, levels, **kwargs)

    return contours
        
  def _plot_req_cont(self):
    """
    plot the normalized standard deviation contour
    """

    if self.normalized:
        refstd = 1.0
    else:
        refstd = self.points[0].s_expd

    t = np.linspace(0, np.pi/2)
    r = np.ones_like(t)
    if self.std_ratios is not None:
        std_ratios = self.std_ratios
    else:
        std_ratios = [refstd]

    for ratio in std_ratios:
        self.polar_ax.plot(t, r*ratio,   '--', color='grey')

    # expected point
    self.polar_ax.text(0.08, refstd, 
            self.title_expected, color='k', 
            horizontalalignment='center',verticalalignment='center',
            fontsize=16)
    self.polar_ax.scatter(0, refstd, c='k', s=75)
    
  def _setup_angle_axis(self):
    """
    set the ticks labels etc for the angle axis
    """
    loc = 'top'
    self.ax.axis[loc].set_axis_direction('bottom')  
    self.ax.axis[loc].toggle(ticklabels=True, label=True)
    self.ax.axis[loc].major_ticklabels.set_axis_direction('top')
    self.ax.axis[loc].label.set_axis_direction('top')        
    self.ax.axis[loc].label.set_text(self.title_polar)

  def _setup_x_axis(self):
    """
    set the ticks labels etc for the x axis
    """
    loc = 'left'
    self.ax.axis[loc].set_axis_direction('bottom') 
    self.ax.axis[loc].label.set_text(self.title_xy)

  def _setup_y_axis(self):
    """
    set the ticks labels etc for the y axis
    """
    loc = 'right'
    self.ax.axis[loc].set_axis_direction('top')   
    self.ax.axis[loc].toggle(ticklabels=True)
    self.ax.axis[loc].major_ticklabels.set_axis_direction('left')      
    self.ax.axis[loc].label.set_text(self.title_xy)  
    
  def _setup_axes(self):
    """
    set the ticks labels etc for the angle x and y axes
    """
    self._setup_angle_axis()
    self._setup_x_axis()
    self._setup_y_axis()

class ZoomInTaylorDiagram(object):
  """
    Given predictions and expected numerical data 
    plot the standard deviation of the differences and correlation between
    expected and predicted in a single-quadrant polar plot, with
    r=stddev and theta=arccos(correlation).

    ZoomInTaylorDiagram can zoom in Taylor diagram (Yi Wang, 07/17/2019)
  """
  def __init__(self, fig=None, max_normed_std=2.5, std_ratios=[1], \
          bias_vmin=None, bias_vmax=None, \
          normalized=True,
          cmap=plt.cm.jet, s=80, title_expected=r'expected'):
    """
    Set up Taylor diagram axes. 
    """   
    self.title_polar    = r'Correlation'
    self.title_xy       = r'Normalized Standard Deviation'    
    self.max_normed_std = max_normed_std
    self.s_min          = 0
    self.std_ratios     = std_ratios
    self.bias_vmin      = bias_vmin
    self.bias_vmax      = bias_vmax
    self.normalized     = normalized
    self.cmap           = cmap
    self.s              = s # marker size
    self.title_expected = title_expected
    
    # Define polar coordinate
    tr = PolarAxes.PolarTransform()
    extreme_finder = angle_helper.ExtremeFinderCycle(20,
                                                     20,
                                                     lon_cycle=360,
                                                     lat_cycle=None,
                                                     lon_minmax=None,
                                                     lat_minmax=(0,
                                                                 np.inf),
                                                     )

    grid_locator1 = angle_helper.LocatorDMS(12)

    tick_formatter1 = angle_helper.FormatterDMS()

    grid_helper = GridHelperCurveLinear(tr,
                                        extreme_finder=extreme_finder,
                                        grid_locator1=grid_locator1,
                                        tick_formatter1=tick_formatter1
                                        )

    # figure
    self.fig = fig                                   
    if self.fig is None:
      self.fig = plt.figure()
    
    # setup axes
    ax = SubplotHost(self.fig, 111, grid_helper=grid_helper)
    self.ax = self.fig.add_subplot(ax)

    # set x and y axis 
    self._setup_axes()
                   
    # attach the polar axes and plot the correlation lines
    self.polar_ax = self.ax.get_aux_axes(tr)
    self.polar_ax.plot([np.pi/2.135, np.pi/2.135], [self.s_min, self.max_normed_std*2], color='grey')
    self.polar_ax.plot([np.pi/2.484, np.pi/2.484], [self.s_min, self.max_normed_std*2], color='grey')
    self.polar_ax.plot([np.pi/3,     np.pi/3    ], [self.s_min, self.max_normed_std*2], color='grey')    
    self.polar_ax.plot([np.pi/3.95,  np.pi/3.95 ], [self.s_min, self.max_normed_std*2], color='grey')    
    self.polar_ax.plot([np.pi/6.95,  np.pi/6.95 ], [self.s_min, self.max_normed_std*2], color='grey')

    # Add norm stddev ratio contours
    self._plot_req_cont(self.std_ratios)

    self.points = []

  def add_prediction(self, expected, predicted, predictor_name='', 
                     plot_pt_id='', plot_marker_id='o'):
    """
    Add a prediction/model to the diagram
    """
    this_point = TaylorDiagramPoint(expected, predicted, 
                                    predictor_name, plot_pt_id, plot_marker_id,
                                    self.normalized)   
    self.points.append(this_point)
    
  def plot(self):
    """
    Place all the loaded points onto the figure 
    """
    rs         = []
    std        = []
    thetas     = []
    biases     = []
    names      = []
    point_tags = []
    marker     = []
    for point in self.points:
      rs.append(point.s_expd)
      std.append(point.s_normd)
      thetas.append(np.arccos(point.corrcoef))
      biases.append(point.bias)
      names.append(point.name)
      point_tags.append(point.point_id)
      marker.append(point.marker_id)

    #Plot the angle values and normalized standard deviations
    # Add the point labels and offset them so they’re readable
    if self.bias_vmin is None:
        self.bias_vmin = min(biases)
    if self.bias_vmax is None:
        self.bias_vmax = max(biases)
    norm = mcolors.Normalize(vmin=self.bias_vmin,vmax=self.bias_vmax)
    for i in range(len(marker)):
        self.sc = self.polar_ax.scatter([thetas[i]], [std[i]], c=[biases[i]], \
                s=self.s, cmap=self.cmap, norm=norm, marker=marker[i], \
                edgecolors='k')
    for i, tag in enumerate(point_tags):                         
      self.polar_ax.text(thetas[i]+0.01, std[i]+0.01, tag,
                         horizontalalignment='center',
                         verticalalignment='bottom', fontsize=16)

  def bias_colorbar(self, **kwargs):
      """
      add bias colorbar
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      cbar = plt.colorbar(self.sc, **kwargs)

      cbar.set_label('Relative bias [%]', fontsize=12)

      return cbar

  def stringID_legend(self, stringID_labels, **kwargs):
      """
      add stringID legend
      stringID_labels is tuple, and every element in the tuple is also a tuple which
      consists marker and label.
      for example, stringID_labels: ( (r'$1$', 'OMI'), (r'$2$', 'OMPS') )
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      handles = []
      for stID_la in stringID_labels:
          stID = stID_la[0]
          la   = stID_la[1]
          handles.append( plt.scatter([], [], marker=stID, label=la, s=self.s, c='k') )

      legend = plt.legend(handles=handles, **kwargs)
      self.ax.add_artist(legend)

  def marker_legend(self, marker_labels, **kwargs):
      """
      add marker legend:
      marker_labels is a tuple, and every element in the tuple is also a tuple which 
      consists marker and label.
      for example, marker_labels: ( ('o', 'prior'), ('s', 'posterior') )
      (Yi Wang, 10/07/2018)
      """

      plt.sca(self.ax)

      handles = []
      for ma_la in marker_labels:
          ma = ma_la[0]
          la = ma_la[1]
          handles.append( plt.scatter([], [], marker=ma, label=la, edgecolors='k', s=self.s, c='w') )

      legend = plt.legend(handles=handles, **kwargs)
      self.ax.add_artist(legend)
                 
  def add_contours(self, levels=5, **kwargs):
    """
    Add constant centered RMS difference contours around reference point
    """

    rs,ts = np.meshgrid(np.linspace(self.s_min,self.max_normed_std*2), \
            np.linspace(0,np.pi/2))

    # Compute centered RMS difference
    if self.normalized:
        refstd = 1.0
    else:
        # not mormailized by expected standard deviation (s_expd)
        # (Yi Wang, 11/02/2019)
        refstd = self.points[0].s_expd
    rms = np.sqrt((refstd)**2 + rs**2 - 2*(refstd)*rs*np.cos(ts))
    
    contours = self.polar_ax.contour(ts, rs, rms, levels, **kwargs)

    return contours
        
  def _plot_req_cont(self, std_ratios):
    """
    plot the normalized standard deviation contour
    """
    self.polar_ax.text(0.08, 1, self.title_expected, color='k', horizontalalignment='center',verticalalignment='center',\
            fontsize=16)
    self.polar_ax.scatter(0, 1, c='k', s=75)

    t = np.linspace(0, np.pi/2)
    r = np.ones_like(t)
    for ratio in std_ratios:
        self.polar_ax.plot(t, r*ratio,   '--', color='grey')  

  def _setup_x_axis(self):
    """
    set the ticks labels etc for the x axis
    """
    self.ax.set_xlim([self.s_min, self.max_normed_std])
    self.ax.set_xlabel(self.title_xy)
    self.ax.axis["bottom"].get_helper().nth_coord_ticks = 1

  def _setup_y_axis(self):
    """
    set the ticks labels etc for the y axis
    """
    self.ax.set_ylim([self.s_min, self.max_normed_std])
    self.ax.set_ylabel(self.title_xy)

    
  def _setup_axes(self):
    """
    set the ticks labels etc for the angle x and y axes
    """
    self.ax.set_aspect('equal')

    self.ax.axis["right"].major_ticklabels.set_visible(False)
    self.ax.axis["top"].major_ticklabels.set_visible(False)

    self._setup_x_axis()
    self._setup_y_axis()



if '__main__'== __name__:

    import copy

    #############################
    # test data
    #############################
    ref = np.array([ 61.37, 38.12, 21.57, 33.25, 32.85, \
                     26.68, 42.07, 44.31, 25.57, 48.34, \
                     53.48, 36.50, 44.94, 60.17, 49.94, \
                     25.67, 39.48, 31.53, 47.71, 33.04, \
                     26.99])

    model1 = np.array([19.47, 21.92,  3.31, 27.10, 34.13, \
                        3.93, 25.44, 24.77,  6.00, 30.17, \
                       34.17, 12.85, 18.76, 29.17, 26.41, \
                        7.06, 23.55, 16.34, 33.52,  9.56, \
                       10.59])

    model2 = np.array([15.68, 18.32,  3.20, 22.11, 19.81, \
                        3.59, 23.96, 25.74,  4.82, 33.47, \
                       30.49, 10.14, 16.21, 31.29, 30.56, \
                        7.32, 13.67, 15.15, 29.49,  8.94, \
                       10.53])

    model1_ratio2 = model1 * 1.2
    model2_ratio2 = model2 * 1.2

    model1_ratio3 = model1 * 2
    model2_ratio3 = model2 * 2

    marker_labels = ( ('o', 'Circle'), \
                      ('s', 'Square') )

    stringID_labels = ( (r'$1$', 'Label 1'), \
                        (r'$2$', 'Label 2'), \
                        (r'$3$', 'Label 3') )

    # zoom in region
    xlim = [0.4, 1.7]
    ylim = [0.4, 1.7]

    # normalized standard deviation line
    std_ratios = [1, 2]
    std = None

    # maximum of normalized standard deviation
    max_normed_std = 2.4
    max_std = 25

    # normlaized bias limit
    bias_vmin = -100
    bias_vmax =  100

    # normalized bias colormap
    cmap = plt.get_cmap('seismic')

    # normalized debiased RMSE levels
    levels = [0.5, 1.0, 1.5, 2.0]
    levels_no_norm = [5, 10, 15, 20]

    ########################################
    # ModTaylorDiagram example (normalized)
    ########################################

    fig1 = plt.figure(figsize=(7,7))
    plt.rcParams.update({'font.size': 12})

    # Initialized an instance
    mtd1 = ModTaylorDiagram(fig=fig1,max_normed_std=max_normed_std, \
            bias_vmin=bias_vmin, bias_vmax=bias_vmax, \
            std_ratios=std_ratios, cmap=cmap, title_expected='expected')

    # add points
    mtd1.add_prediction(ref, model1,        r'', '1', 'o')
    mtd1.add_prediction(ref, model2,        r'', '1', 's')
    mtd1.add_prediction(ref, model1_ratio2, r'', '2', 'o')
    mtd1.add_prediction(ref, model2_ratio2, r'', '2', 's')
    mtd1.add_prediction(ref, model1_ratio3, r'', '3', 'o')
    mtd1.add_prediction(ref, model2_ratio3, r'', '3', 's')
    mtd1.plot()

    # normzliaed debiased RMSE countours
    contours1 = mtd1.add_contours(levels=levels, colors='0.5')
    plt.clabel(contours1, inline=1, fontsize=12, fmt='%1.1f')

    # marker legend
    mtd1.marker_legend( marker_labels )

    # label legend
    mtd1.stringID_legend( stringID_labels, loc='upper left' )

    # plot zoom in region
    mtd1.ax.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], \
                 [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], \
                 'r--', lw=2)

    # normalized bias color bar
    cbar1 = mtd1.bias_colorbar(orientation='horizontal', shrink=0.77, pad=0.08, extend='both')

    ############################################
    # ZoomInTaylorDiagram example normalized
    ############################################

    fig2 = plt.figure(figsize=(7,7))
    plt.rcParams.update({'font.size': 12})

    # Initialized an instance
    mtd2 = ZoomInTaylorDiagram(fig=fig2,max_normed_std=max_normed_std, \
            bias_vmin=bias_vmin, bias_vmax=bias_vmax, \
            std_ratios=std_ratios, cmap=cmap, title_expected='')

    # add points
#    mtd2.add_prediction(ref, model1,        r'', '1', 'o')
#    mtd2.add_prediction(ref, model2,        r'', '1', 's')
#    mtd2.add_prediction(ref, model1_ratio2, r'', '2', 'o')
#    mtd2.add_prediction(ref, model2_ratio2, r'', '2', 's')
#    mtd2.add_prediction(ref, model1_ratio3, r'', '3', 'o')
#    mtd2.add_prediction(ref, model2_ratio3, r'', '3', 's')
    mtd2.points = copy.deepcopy( mtd1.points )
    mtd2.plot()

    # normzliaed debiased RMSE countours
    contours2 = mtd2.add_contours(levels=levels, colors='0.5')
    manual_locations2 = [(1, 0.5), (1, 1.0), (1, 1.5), (1, 2.0)]
    mtd2.polar_ax.clabel(contours2, inline=1, fontsize=12, fmt='%1.1f')

    # marker legend
    mtd2.marker_legend( marker_labels )

    # label legend
    mtd2.stringID_legend( stringID_labels, loc='upper left' )

    # normalized bias color bar
    cbar2 = mtd2.bias_colorbar(orientation='horizontal', shrink=0.77, pad=0.08, extend='both')

    # Zoom in
    mtd2.ax.set_xlim(xlim)
    mtd2.ax.set_ylim(ylim)

    ############################################
    # ModTaylorDiagram example (NOT normalized)
    ############################################

    fig3 = plt.figure(figsize=(7,7))
    plt.rcParams.update({'font.size': 12})

    # Initialized an instance
    mtd3 = ModTaylorDiagram(fig=fig3,max_normed_std=max_std, \
            bias_vmin=bias_vmin, bias_vmax=bias_vmax, \
            normalized=False,
            std_ratios=std, cmap=cmap, title_expected='expected')

    # add points
    mtd3.add_prediction(ref, model1,        r'', '1', 'o')
    mtd3.add_prediction(ref, model2,        r'', '1', 's')
    mtd3.add_prediction(ref, model1_ratio2, r'', '2', 'o')
    mtd3.add_prediction(ref, model2_ratio2, r'', '2', 's')
    mtd3.add_prediction(ref, model1_ratio3, r'', '3', 'o')
    mtd3.add_prediction(ref, model2_ratio3, r'', '3', 's')

    # add standard deviation curve
    mtd3._plot_req_cont()

    mtd3.plot()

    # normzliaed debiased RMSE countours
    contours3 = mtd3.add_contours(levels=levels_no_norm, colors='0.5')
    plt.clabel(contours3, inline=1, fontsize=12, fmt='%2.0f')

    # plot point
    mtd3.plot()

    # marker legend
    mtd3.marker_legend( marker_labels )

    # label legend
    mtd3.stringID_legend( stringID_labels, loc='upper left' )

#    # plot zoom in region
#    mtd3.ax.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], \
#                 [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], \
#                 'r--', lw=2)

    # normalized bias color bar
    cbar3 = mtd3.bias_colorbar(orientation='horizontal', shrink=0.77, pad=0.08, extend='both')

    plt.show()








