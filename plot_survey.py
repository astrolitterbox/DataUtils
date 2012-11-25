# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 19:18:52 2012

@author: astrolitterbox
"""
import numpy as np
from math import sqrt
import itertools
#import db
import matplotlib.pyplot as plt
import matplotlib.font_manager

fig_width_pt = 628.62  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inches
golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in 
fig_height = fig_width*golden_mean+0.3       # height in inches
fig_size = [fig_width,fig_height]    
print fig_size



params = {'backend': 'ps',
          'axes.labelsize': 10,
          'figure.figsize': fig_size,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True, 
          'font': 'serif',
          'font.size': 16,
          'ylabel.fontsize': 20}
plt.rcParams.update(params) 
plt.rc('font', family='serif')
#plt.rc('xtick', labelsize='default')
#plt.rc('ytick', labelsize='default')
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 




class PlotTitles:
    def __init__(self, title, xlabel, ylabel, yMajorTicks, yMinorTicks, xMajorTicks, xMinorTicks):
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.yMajorTicks = yMajorTicks
        self.yMinorTicks = yMinorTicks
        self.xMajorTicks = xMajorTicks
        self.xMinorTicks = xMinorTicks



class GraphData:
    def __init__(self, data, colour, histtype, norm, bins, hatch, legend):
        self.data = data
        self.colour = colour
        self.legends = legend
        self.histtype = histtype
        self.norm = norm
        self.bins = bins
        self.hatch = hatch

                                
class Plots:
    imgDir = 'incl_img/'        
    def plotLogHist(self, graphDataList, filename, plotTitles,  *args):
      
      bins=10**np.arange(*bins)     
      s = plt.figure()
      ax = s.add_subplot(111)
      try:
        args[0]
      except IndexError:
          print 'blah'
      else:  
          v = list(args[0])
          ax.axis(v)
      prop = matplotlib.font_manager.FontProperties(size=10)     
      for gd in graphDataList:
          p1 = plt.hist((gd.data), color=gd.colour, bins = gd.bins, histtype='barstacked', alpha=0.6)    

          plt.legend(gd.legend, loc=gd.loc, markerscale=10, fancybox=True, labelspacing = 0.2, prop=prop, shadow=True)
      plt.title(plotTitles.title)
      plt.xlabel = plotTitles.xlabel
      plt.xscale('log')
      plt.ylabel = plotTitles.ylabel
      
      plt.savefig(self.imgDir+filename)

    def plotHist(self, graphDataList, filename, plotTitles, *args):
      s = plt.figure()

      ax = s.add_subplot(111)
      try:
        args[0]
      except IndexError:
          print 'no axis constraints'
      else:  
          print args[0]
          v = list(args[0])
          ax.axis(v)
      print 'axis', ax.axis
#      if barstacked == True:
#	p1 = plt.hist((gd.data), color=gd.colour, bins = bins, histtype='barstacked', normed=True, alpha=0.6)    
      prop = matplotlib.font_manager.FontProperties(size=16)  
      legendList = []
     
      for gd in graphDataList:
	  legendList.append(gd.legends)
	  print gd.bins
	    
          p1 = plt.hist((gd.data), color=gd.colour, bins = gd.bins, normed=gd.norm, histtype=gd.histtype, edgecolor='black', hatch=gd.hatch, alpha=0.8, linewidth=2)    
	  print gd.legends, 'legend'	          
      plt.legend(legendList, loc=0, markerscale=1, fancybox=False, numpoints = 1, labelspacing = 0.7, prop=prop, shadow=True)
      	
      minor_locator = plt.MultipleLocator(plotTitles.yMinorTicks)
      Ymajor_locator = plt.MultipleLocator(plotTitles.yMajorTicks)  
      major_locator = plt.MultipleLocator(plotTitles.xMajorTicks)      
      Xminor_locator = plt.MultipleLocator(plotTitles.xMinorTicks)   
      ax.xaxis.set_major_locator(major_locator)
      ax.xaxis.set_minor_locator(Xminor_locator)     
      ax.yaxis.set_major_locator(Ymajor_locator)
      ax.yaxis.set_minor_locator(minor_locator)
      #plt.title(plotTitles.title)
      ax.yaxis.label.set_size(20)   	
      ax.xaxis.label.set_size(20) 
      plt.xlabel(plotTitles.xlabel)
      plt.ylabel(plotTitles.ylabel)
      plt.savefig(self.imgDir+filename)
    
    def plotScatter(self, graphDataList, filename, plotTitles, *args):

      s = plt.figure()
      ax = s.add_subplot(111)

      prop = matplotlib.font_manager.FontProperties(size=12)
      legendList = []
      markers = itertools.cycle('.^*')
      for gd in graphDataList:
	  legendList.append(gd.legends)
          p1 = ax.plot(gd.data[0], gd.data[1], linestyle='none', marker=markers.next(), markersize=8, color=gd.colour, mec=gd.colour, alpha = 1)
	  #p1 = ax.plot(gd.data[0], gd.data[1], linewidth=2, markersize=8, color=gd.colour, mec=gd.colour, alpha = 1)	         
	  plt.legend(legendList, loc=2, markerscale=1, numpoints = 1, fancybox=True, labelspacing = 0.8, prop=prop, shadow=True)
      try:
        args[0]
      except IndexError:
          print 'no axis constraints'
      else:
          print 'axis limits', args[0]  
          v = list(args[0])
          ax.axis(v)
      #plt.legend(legendList, loc=0, markerscale=1, fancybox=False, labelspacing = 0.5, prop=prop, shadow=True)
          
      minor_locator = plt.MultipleLocator(plotTitles.yMinorTicks)
      Ymajor_locator = plt.MultipleLocator(plotTitles.yMajorTicks)  
      major_locator = plt.MultipleLocator(plotTitles.xMajorTicks)      
      Xminor_locator = plt.MultipleLocator(plotTitles.xMinorTicks)   
      ax.xaxis.set_major_locator(major_locator)
      ax.xaxis.set_minor_locator(Xminor_locator)     
      ax.yaxis.set_major_locator(Ymajor_locator)
      ax.yaxis.set_minor_locator(minor_locator)
      #plt.title(plotTitles.title)
      ax.yaxis.label.set_size(20)   	
      ax.xaxis.label.set_size(20) 
      plt.xlabel(plotTitles.xlabel)
      plt.ylabel(plotTitles.ylabel)
      plt.savefig(self.imgDir+filename)
