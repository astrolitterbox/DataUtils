# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 19:18:52 2012

@author: astrolitterbox
"""
import numpy as np
#import db
import matplotlib.pyplot as plt
import matplotlib.font_manager
 

class PlotTitles:
    def __init__(self, title, xlabel, ylabel):
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel



class GraphData:
    def __init__(self, data, colour, legend):
        self.data = data
        self.colour = colour
        self.legends = legend

class Plots:
    imgDir = './img/'        
    def plotLogHist(self, graphDataList, filename, plotTitles, bins, *args):
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
      prop = matplotlib.font_manager.FontProperties(size=8)     
      for gd in graphDataList:
          p1 = plt.hist((gd.data), color=gd.colour, bins = bins, histtype='barstacked', alpha=0.6)    
          print gd.legend
          plt.legend(gd.legend, loc=0, markerscale=10, fancybox=True, labelspacing = 0.2, prop=prop, shadow=True)
      plt.title(plotTitles.title)
      plt.xlabel = plotTitles.xlabel
      plt.xscale('log')
      plt.ylabel = plotTitles.ylabel
      plt.savefig(self.imgDir+filename)

    def plotHist(self, graphDataList, filename, plotTitles, bins, *args):
      s = plt.figure()
      ax = s.add_subplot(111)
      try:
        args[0]
      except IndexError:
          print 'no axis constraints'
      else:  
          v = list(args[0])
          ax.axis(v)
#      if barstacked == True:
#	p1 = plt.hist((gd.data), color=gd.colour, bins = bins, histtype='barstacked', normed=True, alpha=0.6)    
      prop = matplotlib.font_manager.FontProperties(size=8)     
      for gd in graphDataList:
	  
          p1 = plt.hist((gd.data), color=gd.colour, bins = bins, normed=True, histtype='step', alpha=0.75, linewidth=3)    
          print gd.legends            
          plt.legend(gd.legends, loc=0, markerscale=10, fancybox=True, labelspacing = 0.2, prop=prop, shadow=True)
      plt.title(plotTitles.title)
      plt.xlabel = plotTitles.xlabel
      plt.ylabel = plotTitles.ylabel
      plt.savefig(self.imgDir+filename)
    
    def plotScatter(self, graphDataList, filename, plotTitles, *args):
      s = plt.figure()
      ax = s.add_subplot(111)
      try:
        args[0]
      except IndexError:
          print 'no axis constraints'
      else:  
          v = list(args[0])
          ax.axis(v)
      prop = matplotlib.font_manager.FontProperties(size=8)     
      for gd in graphDataList:
          p1 = ax.plot(gd.data[0], gd.data[1], '.', markersize=6, color=gd.colour, mec=gd.colour, alpha = 0.9) 
          #plt.legend([p1[0]], gd.legend,  loc=0, markerscale=10, fancybox=True, labelspacing = 0.2, prop=prop, shadow=True)
      plt.title(plotTitles.title)
      plt.xlabel = plotTitles.xlabel
      plt.ylabel = plotTitles.ylabel
      plt.savefig(self.imgDir+filename)
