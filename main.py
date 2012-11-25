# -*- coding: utf-8 -*-
import numpy as np
import math as math
import matplotlib.pylab as plt
import kcorrect as KC
from kcorrect.sdss import SDSSKCorrect, SDSSPhotoZ, SDSS
from kcorrect.utils.cosmology import ztodm
from scipy import integrate
from scipy.stats import ks_2samp as ks_2samp
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import pyfits as pyfits

#import sami_db
import db as db
import plot_survey as plot
#from sami_new import getDataArray as fullDataArray
#from sami_new import plotScatter as OldSchoolPlot
import matplotlib.font_manager 
import csv
import readAtlas
import collections
from math import sqrt

c = 299792.458 
pi = 3.14159265
#     Cosmological parameters
H0 = 70.0 #km/s/Mpc, 1 Mpc= 3.08568025*1e19 km
Wm = 0.3
Wl = 0.7
Wh = c / H0 
Wk = 1 - Wm - Wl 
tH = (3.08568025*1e19/H0)/(3600 * 365 * 24*10e9) #Hubble time in Gyr
skyArea = 41252.9612494193 #degrees      
dH = c/H0 #in Mpc
imgDir = 'img/'
dbFile = 'CALIFA.sqlite'
dataFile = 'data/Califa.csv'
morphDataFile = 'morph.dat'
observedList = 'list_observed.txt'
db_dataFile = 'db_data.txt'

#colour settings

ASColour = '#666666'
RASColour = '#ED9121'
RDASColour = '#EE0000'
DASColour = '#99FF33'
DRASColour = '#006633'
CALIFAColour = '#1B3F8B'
AIPBlue = '#00417A'

#utilities
def E_z(z): #time derivative of log(a(t)), used for integration
	return 1/(np.sqrt((Wm*(1+z)**3 + Wk*(1+z)**2 + Wl)))

def comovingDistance(z):
	Dc = dH * 1000* integrate.quad(E_z, 0, z)[0] #in kpc
	return Dc

def angular2physical(reff, z): #return physical effective diameter of the galaxy in kpc
    return (math.radians(2*reff/3600) *(comovingDistance(z)) / (1 + z)**2)

def getAbsMag(z, mag, ext):
	print z
	d = comovingDistance(z)
        dm = ztodm(z, (0.3, 0.7, 0.7))
	absmag = mag - dm #- 3.1 * ext
	#print dm, 'dm', absmag, 'absmag'
	return absmag
	
def convert(data):
     tempDATA = []
     for i in data:
     	 #print i
         tempDATA.append([float(j) for j in i])
     return np.asarray(tempDATA)

def findOverlap(a, b):
   a_multiset = collections.Counter(list(a))
   b_multiset = collections.Counter(list(b))
   return list((a_multiset & b_multiset).elements())
  
def getDataArray(data, zData):
  #objid,ra,dec,b,l,petroMag_u,petroMag_g,petroMag_r,petroMag_i,petroMag_z,petroMagErr_u,petroMagErr_g,petroMagErr_r,petroMagErr_i,petroMagErr_z,extinction_u,extinction_g,extinction_r,extinction_i,extinction_z,isoA_r,isoB_r,z,rho,lnLExp_r,expRad_r,lnLDeV_r,deVRad_r
  dataArray = np.zeros((len(data), 27))
  dataArray[:, 0] = data[0]['objid']
  dataArray[:, 1] = data[:]['ra']
  dataArray[:, 2] = data[:]['dec']
  dataArray[:, 3] = data[:]['petroMag_u']
  dataArray[:, 4] = data[:]['petroMag_g']
  dataArray[:, 5] = data[:]['petroMag_r']
  dataArray[:, 6] = data[:]['petroMag_i']
  dataArray[:, 7] = data[:]['petroMag_z']
  dataArray[:, 8] = data[:]['petroMagErr_u']
  dataArray[:, 9] = data[:]['petroMagErr_g']
  dataArray[:, 10] = data[:]['petroMagErr_r']
  dataArray[:, 11] = data[:]['petroMagErr_i']
  dataArray[:, 12] = data[:]['petroMagErr_z']
  dataArray[:,13] = data[:]['extinction_u']
  dataArray[:,14] = data[:]['extinction_g']
  dataArray[:,15] = data[:]['extinction_r']
  dataArray[:,16] = data[:]['extinction_i']
  dataArray[:,17] = data[:]['extinction_z']
  dataArray[:,18] = data[:]['isoA_r']
  dataArray[:,19] = data[:]['isoB_r']
  dataArray[:,20] = zData
  
  
  redshifts = dataArray[:,20]
  maggies = dataArray[:,3:8]
  maggies_err = dataArray[:,8:13]
  extinction = dataArray[:,13:18]
  kc = SDSSKCorrect(redshifts, maggies, maggies_err, extinction, cosmo=(Wm, Wl, H0/100))
  kcorr = kc.kcorrect()
  absmag = kc.absmag()
  dataArray[:,21] = absmag[:, 2]
  dataArray[:,22] = data[:]['petroMag_g'] - data[:]['extinction_g'] - (data[:]['petroMag_r'] - data[:]['extinction_r'])
  sb = data[:]['rho'] + data[:]['petroMag_r']
  dataArray[:,23] = sb
  stMass = np.zeros((len(data), 1))
  size = np.zeros((len(data), 1))
  reff = np.zeros((len(data), 1))
  for i in range(0, len(data)):
    if data[i]['lnLDeV_r'] >= data[i]['lnLExp_r']:
      reff[i] = data[i]['deVRad_r']
    else:
      reff[i] = data[i]['expRad_r']
    size[i] = angular2physical(reff[i], redshifts[i])
    dataArray[i,24] = reff[i]
    dataArray[i,25] = size[i]
  coeffs = kc.coeffs
  tmremain = np.array([[0.601525, 0.941511, 0.607033, 0.523732, 0.763937]])
  ones = np.ones((1, len(redshifts)))
  prod = np.dot(tmremain.T, ones).T 
  modelMasses = coeffs*prod
  #print modelMasses.shape
  mass = np.sum(modelMasses, axis=1)
  for i in range (0, (len(data))):
    distmod = KC.utils.cosmology.ztodm(redshifts[i])
    exp = 10 ** (0.4 * distmod)
    stMass[i] = mass[i] * exp
    dataArray[i,26] = stMass[i]
  return dataArray

def createFakeInclinations():
  randomInclinations = (180*np.random.random_sample((2000, 1)))
  
  #print randomInclinations
  randomAxisRatio = np.cos(randomInclinations)
  validAxisRatio = randomAxisRatio[np.where(randomAxisRatio > 0.09)]
  print len(validAxisRatio), min(validAxisRatio)
  return correctForThickness((validAxisRatio))


def correctForThickness(axisRatio):
  print 'correcting', min(axisRatio), 'min'
  return np.degrees(np.arccos(np.sqrt((axisRatio**2 - 0.09**2)/(1 - 0.09**2))))

def list_to_int(l):
	return int("".join(str(x) for x in l))


def main():
  #comparison of Nadine's and SDSS PA, _pixel_coordinates!
  #ic = convert(db.dbUtils.getFromDB('isoPhi_r', 'CALIFA.sqlite', 'mothersample', ' where CALIFA_ID = 544'))
  #print ic
  #ic = convert(db.dbUtils.getFromDB('pa_align', 'CALIFA.sqlite', 'nadine', ' where CALIFA_ID = 544'))
  #print ic+180
  #plotting axis ratio distribution for MS, DAS and RAS

  #different diameter change corrections, from Patureli et al.
  '''
  D = 10
  q = 0.2
  incl = np.arange(0, 90)
  axisRatio = np.sqrt(np.cos(np.radians(incl))**2 * (1 - q**2)+q**2)
  C = 0.01
  Dincl_1 = 10**(math.log10(D) + C*np.log10(1/axisRatio))
  Dincl_1 = ((Dincl_1 - D)/(0.01 * D)) #percentage change
  C = 0.05
  Dincl_2 = 10**(math.log10(D) + C*np.log10(1/axisRatio))
  Dincl_2 = ((Dincl_2 - D)/(0.01 * D)) #percentage change
  C = 0.1
  Dincl_3 = 10**(math.log10(D) + C*np.log10(1/axisRatio))
  Dincl_3 = ((Dincl_3 - D)/(0.01 * D)) #percentage change
  graph11 = plot.Plots()
  
  Data1 = plot.GraphData((incl, Dincl_1), AIPBlue, 'step', False, 20, ",", 'C = 0.01')
  Data2 = plot.GraphData((incl, Dincl_2), 'grey', 'step', False, 20, "x", 'C = 0.1')
  Data3 = plot.GraphData((incl, Dincl_3), 'r', 'step', False, 20, ".", 'C = 0.2')

  bins = graph11.plotScatter([Data1, Data2, Data3], 'incl_corr.pdf', plot.PlotTitles('Projected isophotal diameter change with inclination',  'Inclination, deg', 'Projected $D_{iso}$ change, $\%$', 5, 2, 10, 5))
  '''

 #angular diameter isoA_r vs. app.mag. distribution
  isoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  z = convert(db.dbUtils.getFromDB('z',  'CALIFA.sqlite', 'mothersample'))
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  
  Incl = convert(db.dbUtils.getFromDB('isoB_r',  'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample'))
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  #data for SDSS
  
  absmag_fo = r[faceOn]
  absmag_mi = r[midincl]
  absmag_eo = r[edgeOn]
  
  size_fo = isoA_r[faceOn]
  size_mi = isoA_r[midincl]
  size_eo = isoA_r[edgeOn]
 
  graph11 = plot.Plots()
  fo_Data = plot.GraphData((absmag_fo, size_fo), AIPBlue, 'step', False, 20, ",",  'b/a $>$ 0.7')
  mi_Data = plot.GraphData((absmag_mi, size_mi), 'grey', 'step', False, 20, "x", '0.7 $>$ b/a $>$ 0.3')
  eo_Data = plot.GraphData((absmag_eo, size_eo), 'r', 'step', False, 20, ".", 'b/a $<$ 0.3')
  bins = graph11.plotScatter([mi_Data, fo_Data, eo_Data], 'app_mag_isoA', plot.PlotTitles('Angular size vs. apparent r band magnitude',  'SDSS Petrosian r magnitude, mag', 'Angular isophotal radius, arcsec', 20, 5, 2, 0.5), (11, 16, 42, 80))




  
  #physical size from isoA_r vs. absmag distribution
  isoA_r = convert(db.dbUtils.getFromDB('isoA', 'CALIFA.sqlite', 'lucie_ellipse'))
  print isoA_r.shape
  z = convert(db.dbUtils.getFromDB('z',  'CALIFA.sqlite', 'mothersample'))[0:937]
  #r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  r = convert(db.dbUtils.getFromDB('el_mag', 'CALIFA.sqlite', 'gc'))[0:937]
  print r.shape

  size = np.empty((len(z), 1))
  for i in range(0, len(z)):
  	  size[i] = angular2physical(isoA_r[i], z[i])
  
  
  absmag = np.empty((len(z), 1))
  for i in range(0, len(z)):

  	  absmag[i] = getAbsMag(z[i], r[i], 0)
  Incl = convert(db.dbUtils.getFromDB('ba',  'CALIFA.sqlite', 'nadine'))[0:937]
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  #data for SDSS
  
  absmag_fo = absmag[faceOn]
  absmag_mi = absmag[midincl]
  absmag_eo = absmag[edgeOn]
  
  size_fo = size[faceOn]
  size_mi = size[midincl]
  size_eo = size[edgeOn]
 
  graph11 = plot.Plots()
  fo_Data = plot.GraphData((absmag_fo, size_fo), AIPBlue, 'step', False, 20, ",",  'b/a $>$ 0.7')
  mi_Data = plot.GraphData((absmag_mi, size_mi), 'grey', 'step', False, 20, "x", '0.7 $>$ b/a $>$ 0.3')
  eo_Data = plot.GraphData((absmag_eo, size_eo), 'r', 'step', False, 20, ".", 'b/a $<$ 0.3')
  bins = graph11.plotScatter([mi_Data, fo_Data, eo_Data], 'sl_rel_gc_vs_lucie', plot.PlotTitles('Size vs. absolute magnitude',  'Growth curve Mr, mag', 'Linear isophotal size from ellipse fit, kpc', 20, 5, 2, 0.5), (-16, -24, 5, 80))
  exit()
  #SDSS and Lucie's b/a comparison
  from matplotlib  import cm
  dbDir = ""
  graph = plot.Plots()
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))[0:937]/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))[0:937]

  lucie_ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'lucie'))  #parsing tuples 
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Lucie isophotal $b/a$ vs. SDSS b/a ",fontsize=14)
  ax.set_xlabel("SDSS b/a",fontsize=12)
  ax.set_ylabel("Isophotal b/a",fontsize=12)
  sc = ax.scatter(isoA_r_axis_ratio,lucie_ba,s=6,c='black', marker = 'o', edgecolor='None')
  ax.axis([0, 1, 0, 1])
  #cb = fig.colorbar(sc)
  #cb.set_label('b/a')
  plt.savefig('incl_img/b_a_comparison')
  #plotData = plot.GraphData(((sdss_hlr, lucie_hlr)), Incl, 'step', False, 20, ".", 'SDSS $r_{50}$ vs. Lucie isophotal $r_e$')    
  #graph.plotScatter([plotData], "circ_hlr_vs_sdss", plot.PlotTitles("Comparison between Lucie's and SDSS HLR values", "SDSS Petrosian R50, arcsec", "Lucies isophotal $r_e$, arcsec", 20, 5, 20, 5), (0, 50, 0, 45))
  
  
  
  
  #SDSS and Lucie's comparison
  from matplotlib  import cm
  dbDir = ""
  graph = plot.Plots()
  sdss_hlr = convert(db.dbUtils.getFromDB('petroR50_r', dbDir+'CALIFA.sqlite', 'mothersample'))[0:937]
  Incl = (convert(db.dbUtils.getFromDB('isoB_r',  'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample')))[0:937]
  lucie_hlr = convert(db.dbUtils.getFromDB('re', dbDir+'CALIFA.sqlite', 'lucie'))  #parsing tuples 
  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Lucie isophotal $r_e$ vs. SDSS $PetroR_{50}$ ",fontsize=14)
  ax.set_xlabel("SDSS Petrosian half-light radius, arcsec",fontsize=12)
  ax.set_ylabel("Isophotal half-light radius, arcsec",fontsize=12)
  sc = ax.scatter(sdss_hlr,lucie_hlr,s=6,c=Incl, marker = 'o', edgecolor='None', cmap = cm.jet)
  ax.axis([0, 50, 0, 50])
  cb = fig.colorbar(sc)
  cb.set_label('b/a')
  plt.savefig('incl_img/lucie_vs_sdss')
  #plotData = plot.GraphData(((sdss_hlr, lucie_hlr)), Incl, 'step', False, 20, ".", 'SDSS $r_{50}$ vs. Lucie isophotal $r_e$')    
  #graph.plotScatter([plotData], "circ_hlr_vs_sdss", plot.PlotTitles("Comparison between Lucie's and SDSS HLR values", "SDSS Petrosian R50, arcsec", "Lucies isophotal $r_e$, arcsec", 20, 5, 20, 5), (0, 50, 0, 45))

  #axis ratio distribution
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'Spirals'))

  graph = plot.Plots()
  hist, bins = np.histogram(spiralsIncl, bins = 10)

  width=1*(bins[1]-bins[0])
  #center=(bins[:-1]+bins[1:])/2

  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio), 'none', 'bar', False, 1/width, "o",  'Full CALIFA sample')  
  spiralsIncl = plot.GraphData((spiralsIncl), 'red', 'stepfilled', False, 1/width, "//",  'Spirals')
  graph.plotHist([spiralsIncl, isoA_r_AxisRatioData], 'ba_dist.pdf', plot.PlotTitles("b/a distribution of the whole MS and of the spiral galaxies", "b/a", "N",  20, 10, 0.2, 0.1)) 
  
  #ancestor and CALIFA MS r mag vs. b/a:
  r_DAS = convert(db.dbUtils.getFromDB('petroMag_r', 'ancestor_samples.sqlite', 'DAS', ' where petroR90_r/petroR50_r < 2.5'))
  r_RAS = convert(db.dbUtils.getFromDB('petroMag_r', 'ancestor_samples.sqlite', 'RAS', ' where isoA_r > 20 and petroR90_r/petroR50_r < 2.5'))
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  isoA_r_DAS = convert(db.dbUtils.getFromDB('isoB_r', 'ancestor_samples.sqlite', 'DAS', ' where petroR90_r/petroR50_r < 2.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'ancestor_samples.sqlite', 'DAS', ' where petroR90_r/petroR50_r < 2.5'))
  isoA_r_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'ancestor_samples.sqlite', 'RAS', ' where isoA_r > 20 and petroR90_r/petroR50_r < 2.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'ancestor_samples.sqlite', 'RAS', ' where isoA_r > 20 and petroR90_r/petroR50_r < 2.5'))

  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  
  graph = plot.Plots()
  hist, bins = np.histogram(isoA_r_axis_ratio, bins = 10)

  width=1*(bins[1]-bins[0])
  #center=(bins[:-1]+bins[1:])/2

  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio, r), 'k', 'bar', True, 1/width, "o",  'Petrosian r, CALIFA sample')  
  DAS_AxisRatioData = plot.GraphData((isoA_r_DAS, r_DAS), 'grey', 'stepfilled', True, 1/width, "", 'Petrosian r, diameter-selected')  
  RAS_AxisRatioData = plot.GraphData((isoA_r_RAS, r_RAS), 'red', 'stepfilled', True, 1/width, "//", 'Petrosian r, redshift-selected')

  graph.plotScatter([RAS_AxisRatioData, DAS_AxisRatioData, isoA_r_AxisRatioData], '3_samples_mag_vs_ba', plot.PlotTitles("r band magnitude vs. b/afor the mothersample, DAS and RAS", "b/a", "r band SDSS magnitude", 1, 0.5, 0.2, 0.1))  
  
  
  
  
  exit()
  #3 samples
  isoA_r_DAS = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DAS', ' where petroR90_r/petroR50_r < 2.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS', ' where petroR90_r/petroR50_r < 2.5'))
  isoA_r_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS2', ' where isoA_r > 20 and petroR90_r/petroR50_r < 2.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS2', ' where isoA_r > 20 and petroR90_r/petroR50_r < 2.5'))
  print isoA_r_RAS.shape 
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  
  graph = plot.Plots()
  hist, bins = np.histogram(isoA_r_axis_ratio, bins = 10)

  width=1*(bins[1]-bins[0])
  #center=(bins[:-1]+bins[1:])/2

  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio), 'white', 'bar', True, 1/width, "o",  'b/a, CALIFA sample')  
  DAS_AxisRatioData = plot.GraphData((isoA_r_DAS), 'grey', 'stepfilled', True, 1/width, "", 'b/a, diameter-selected')  
  RAS_AxisRatioData = plot.GraphData((isoA_r_RAS), 'red', 'stepfilled', True, 1/width, "//", 'b/a, redshift-selected')

  graph.plotHist([RAS_AxisRatioData, DAS_AxisRatioData, isoA_r_AxisRatioData], '3_samples.pdf', plot.PlotTitles("b/a distributions for the mothersample, DAS and RAS", "b/a", "N, normalised number", 0.5, 0.1, 0.2, 0.1))  
  
    
 
  
  '''
  stMass = convert(db.dbUtils.getFromDB('StMass', 'CALIFA.sqlite', 'mothersample'))
  graph11 = plot.Plots()
  plot_Data = plot.GraphData((stMass), 'k', 'step', False, 20, 'Galaktikų žvaigždžių masė')
  bins = graph11.plotLogHist([plot_Data], 'st_masses', plot.PlotTitles("Galaktikų žvaigždžių masių pasiskirstymas",  'M, M_sol', ''))
  exit()
  err = convert(db.dbUtils.getFromDB('PA_err', 'CALIFA.sqlite', 'nadine')) 
  goodPoints = np.where(err < 20)
  print err, err.shape, 'err'
  pa = convert(db.dbUtils.getFromDB('pa_align', 'CALIFA.sqlite', 'nadine'))
  pa_sdss = convert(db.dbUtils.getFromDB('isoPhi_r', 'CALIFA.sqlite', 'mothersample'))[0:937]
  pa[np.where(pa < 0)] = pa[np.where(pa < 0)] + 180
  
  print pa.shape, pa_sdss.shape, 'length of pa vectors'
  
  Incl = (convert(db.dbUtils.getFromDB('isoB_r',  'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample')))[0:937]
  print Incl.shape
  Incl_Nadine = convert(db.dbUtils.getFromDB('ba',  'CALIFA.sqlite', 'nadine'))
 
  graph11 = plot.Plots()
  plot_Data = plot.GraphData((Incl, pa-pa_sdss), 'k', 'step', False, 20, 'PA_gc - PA_sdss')
  
  bins = graph11.plotScatter([plot_Data], 'pa_diff_vs_sdss_ba_all', plot.PlotTitles("PA difference vs. axis ratio",  'Axis ratio, SDSS isophotal', 'Difference between PA values, deg'))
  
  
  plot_Data = plot.GraphData((Incl_Nadine, pa-pa_sdss), 'k', 'step', False, 20, 'PA_gc - PA_sdss')
  bins = graph11.plotScatter([plot_Data], 'pa_diff_vs_growth_curve_ba_all', plot.PlotTitles("PA difference vs. axis ratio",  'Axis ratio, GC analysis', 'Difference between PA values, deg'), (0, 1, -200, 200))
  print Incl.shape, 'incl_shape'
  pa_sdss = pa_sdss[goodPoints]
  pa = pa[goodPoints]
  Incl = Incl[goodPoints]
  Incl_Nadine = Incl_Nadine[goodPoints]
  print Incl_Nadine.shape, 'Nincl_shape'
  #b/a for SDSS
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  #data for SDSS
  
  faceonData = pa[faceOn]
  midinclData = pa[midincl]
  edgeOnData = pa[edgeOn]
  
  
  
  faceonsdss = pa_sdss[faceOn]
  midinclsdss = pa_sdss[midincl]
  edgeOnsdss = pa_sdss[edgeOn]

  graph11 = plot.Plots()
  FO_Data = plot.GraphData((faceonsdss, faceonData), 'b', 'step', False, 20, 'b/a > 0.7')
  MI_Data = plot.GraphData((midinclsdss, midinclData), 'k', 'step', False, 20, '0.7 > b/a >0.3')
  EO_Data = plot.GraphData((edgeOnsdss, edgeOnData), 'r', 'step', False, 20, 'b/a < 0.3')
  
  plot_Data = plot.GraphData((Incl, pa-pa_sdss), 'k', 'step', False, 20, 'PA_gc - PA_sdss')
  
  bins = graph11.plotScatter([plot_Data], 'pa_diff_vs_sdss_ba', plot.PlotTitles("PA difference vs. axis ratio",  'Axis ratio, SDSS isophotal', 'Difference between PA values, deg'))
  
  bins = graph11.plotScatter([FO_Data, MI_Data, EO_Data], 'pa_growth_vs_sdss', plot.PlotTitles("Comparison between the SDSS and growth curve PA values",  'PA, SDSS, deg', 'PA, by Nadine, deg'))
  
  #for Nadines ba:
  
  
  plot_Data = plot.GraphData((Incl_Nadine, pa-pa_sdss), 'k', 'step', False, 20, 'PA_gc - PA_sdss')
  
  bins = graph11.plotScatter([plot_Data], 'pa_diff_vs_growth_curve_ba', plot.PlotTitles("PA difference vs. axis ratio",  'Axis ratio, GC analysis', 'Difference between PA values, deg'), (0, 1, -200, 200))
  
  
  

  #b/a for Nadine
  faceOn = np.where(Incl_Nadine > 0.7)
  midincl = np.where((Incl_Nadine > 0.3) & (Incl_Nadine < 0.7))
  edgeOn = np.where((Incl_Nadine < 0.3))
  
  faceonData = pa[faceOn]
  midinclData = pa[midincl]
  edgeOnData = pa[edgeOn]
  
  faceonsdss = pa_sdss[faceOn]
  midinclsdss = pa_sdss[midincl]
  edgeOnsdss = pa_sdss[edgeOn]

  graph11 = plot.Plots()
  FO_Data = plot.GraphData((faceonsdss, faceonData), 'b', 'step', False, 20, 'b/a > 0.7')
  MI_Data = plot.GraphData((midinclsdss, midinclData), 'k', 'step', False, 20, '0.7 > b/a >0.3')
  EO_Data = plot.GraphData((edgeOnsdss, edgeOnData), 'r', 'step', False, 20, 'b/a < 0.3')
  

  bins = graph11.plotScatter([FO_Data, MI_Data, EO_Data], 'pa_growth_vs_sdss_nadine_ba', plot.PlotTitles("Comparison between the SDSS and growth curve PA values",  'PA, SDSS, deg', 'PA, by Nadine, deg'))
  exit()
  #data output for CMD
  absmag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'mothersample'))
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  g = convert(db.dbUtils.getFromDB('petroMag_g', 'mothersample', 'mothersample'))
  gr = g - r
  
  Oabsmag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'observed'))
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'observed'))
  g = convert(db.dbUtils.getFromDB('petroMag_g', 'mothersample', 'observed'))
  Ogr = g - r
  data = np.array((gr, absmag))
  np.savetxt("Mr.txt", absmag)
  np.savetxt("gr.txt", gr)

  np.savetxt("OMr.txt", Oabsmag)
  np.savetxt("Ogr.txt", Ogr)


  exit()
  #isoA_r comparison for original and mothersample tables
  isoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  isoA_rOrig = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'Final'))
  print np.where(isoA_r < 45), np.where(isoA_rOrig < 45)
  print min(isoA_r)[0:937], min(isoA_rOrig)
  print isoA_rOrig[303]
  graph11 = plot.Plots()
  isoA_Data = plot.GraphData((isoA_r[0:937], isoA_rOrig), 'k', 'step', False, 20, 'IsoA_r')
  bins = graph11.plotScatter([isoA_Data], 'isoA_MS_vs_original', plot.PlotTitles("Comparison between the original and derived MS table",  'isoA_r, arcsec', 'isoA_r, arcsec'))


  #m_r comparison for original and mothersample tables
  m_r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  m_rOrig = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'Final'))
  graph11 = plot.Plots()
  mr_Data = plot.GraphData((m_r[0:937], m_rOrig), 'k', 'step', False, 20, 'm_r')
  bins = graph11.plotScatter([mr_Data], 'm_r_MS_vs_original', plot.PlotTitles("Comparison between the original and derived m_r",  'm_r, mag', 'm_r, mag'))

  exit()
  objidOrig = db.dbUtils.getFromDB('ra, dec', 'CALIFA.sqlite', 'original')[0:10]
  objidMS = db.dbUtils.getFromDB('ra, dec', 'CALIFA.sqlite', 'mothersample')[0:10]
  print objidOrig, objidMS
  exit()
  #physical size from isoA_r vs. absmag distribution, compared with RAS
  isoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  z = convert(db.dbUtils.getFromDB('z',  'CALIFA.sqlite', 'mothersample'))
  #ext = convert(db.dbUtils.getFromDB('extinction_r',  'CALIFA.sqlite', 'mothersample'))
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  size = np.empty((len(z), 1))
  for i in range(0, len(z)):
          print i
  	  size[i] = angular2physical(isoA_r[i], z[i])
  
  absmag = np.empty((len(z), 1))
  for i in range(0, len(z)):
  	  absmag[i] = getAbsMag(z[i], r[i], 0)
  
  Incl = convert(db.dbUtils.getFromDB('isoB_r',  'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample')) 
  
  RAS_Incl = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'random'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'random')) 
  
  RAS_r = convert(db.dbUtils.getFromDB('petroMag_r', 'RAS', 'random'))
  RAS_ext = convert(db.dbUtils.getFromDB('extinction_r', 'RAS', 'random'))
  RAS_z = convert(db.dbUtils.getFromDB('z', 'RAS', 'random'))
  RAS_isoA_r = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'random'))

  RAS_size = np.empty((len(RAS_z), 1))
  for i in range(0, len(RAS_z)):
          print i
  	  RAS_size[i] = angular2physical(RAS_isoA_r[i], RAS_z[i])
  
  RAS_Mr = np.empty((len(RAS_z), 1)) 
  for i in range(0, len(RAS_z)):
  	print i
	RAS_Mr[i]  = getAbsMag(RAS_z[i], (RAS_r[i]), RAS_ext[i])
  
  #b/a for CALIFA
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnSize = size[faceOn]
  midInclsSize = size[midincl]
  edgeOnSize = size[edgeOn]


  faceOnMr = absmag[faceOn]
  midInclsMr = absmag[midincl]
  edgeOnInclsMr = absmag[edgeOn]
  
  #b/a for DAS
  faceOn = np.where(RAS_Incl > 0.7)
  midincl = np.where((RAS_Incl > 0.3) & (RAS_Incl < 0.7))
  edgeOn = np.where((RAS_Incl < 0.3))
  
  faceOnDASSize = RAS_size[faceOn]
  midInclsDASSize = RAS_size[midincl]
  edgeOnDASSize = RAS_size[edgeOn]


  faceOnDASMr = RAS_Mr[faceOn]
  midInclsDASMr = RAS_Mr[midincl]
  edgeOnInclsDASMr = RAS_Mr[edgeOn]
  
  
  graph11 = plot.Plots()
  FOsize_Data = plot.GraphData((faceOnMr, faceOnSize), 'b', 'step', False, 20, 'CALIFA b/a > 0.7')
  MIsize_Data = plot.GraphData((midInclsMr, midInclsSize), 'k', 'step', False, 20, 'CALIFA b/a between 0.3 and 0.7')
  EOsize_Data = plot.GraphData((edgeOnInclsMr, edgeOnSize), 'r', 'step', False, 20, 'CALIFA b/a < 0.3')

  DAS_FOsize_Data = plot.GraphData((faceOnDASMr, faceOnDASSize), DASColour, 'step', False, 20, 'RAS b/a > 0.7')
  DAS_MIsize_Data = plot.GraphData((midInclsDASMr, midInclsDASSize), DRASColour, 'step', False, 20, 'RAS b/a between 0.3 and 0.7')
  DAS_EOsize_Data = plot.GraphData((edgeOnInclsDASMr, edgeOnDASSize), 'grey', 'step', False, 20, 'RAS b/a < 0.3')


  bins = graph11.plotScatter([DAS_FOsize_Data, DAS_MIsize_Data, DAS_EOsize_Data, FOsize_Data, MIsize_Data, EOsize_Data], 'iso_size_vs_absmag_RAS', plot.PlotTitles("Physical size vs. absolute magnitude",  'M_r, mag', 'Isophotal size, kpc'), (-16, -24, 0, 90))



  #physical size from isoA_r vs. absmag distribution, compared with DAS
  
  
  isoA_r = convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample'))
  z = convert(db.dbUtils.getFromDB('z',  'CALIFA.sqlite', 'mothersample'))
  ext = 0#convert(db.dbUtils.getFromDB('extinction_r',  'CALIFA.sqlite', 'mothersample'))
  r = convert(db.dbUtils.getFromDB('petroMag_r',  'CALIFA.sqlite', 'mothersample'))
  size = np.empty((len(z), 1))
  for i in range(0, len(z)):
          print i
  	  size[i] = angular2physical(isoA_r[i], z[i])
  
  absmag = np.empty((len(z), 1))
  for i in range(0, len(z)):
  	  absmag[i] = getAbsMag(z[i], r[i], 0)
  
  Incl = convert(db.dbUtils.getFromDB('isoB_r',  'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r',  'CALIFA.sqlite', 'mothersample')) 
  
  DAS_Incl = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DAS_z'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS_z')) 
  
  DAS_r = convert(db.dbUtils.getFromDB('petroMag_r', 'DAS', 'DAS_z'))
  DAS_ext = convert(db.dbUtils.getFromDB('extinction_r', 'DAS', 'DAS_z'))
  DAS_z = convert(db.dbUtils.getFromDB('z', 'DAS', 'DAS_z'))
  DAS_isoA_r = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS_z'))

  DAS_size = np.empty((len(DAS_z), 1))
  for i in range(0, len(DAS_z)):
          print i
  	  DAS_size[i] = angular2physical(DAS_isoA_r[i], DAS_z[i])
  
  DAS_Mr = np.empty((len(DAS_z), 1)) 
  for i in range(0, len(DAS_z)):
	DAS_Mr[i]  = getAbsMag(DAS_z[i], (DAS_r[i]), DAS_ext[i])
  
  #b/a for CALIFA
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnSize = size[faceOn]
  midInclsSize = size[midincl]
  edgeOnSize = size[edgeOn]


  faceOnMr = absmag[faceOn]
  midInclsMr = absmag[midincl]
  edgeOnInclsMr = absmag[edgeOn]
  
  #b/a for DAS
  faceOn = np.where(DAS_Incl > 0.7)
  midincl = np.where((DAS_Incl > 0.3) & (DAS_Incl < 0.7))
  edgeOn = np.where((DAS_Incl < 0.3))
  
  faceOnDASSize = DAS_size[faceOn]
  midInclsDASSize = DAS_size[midincl]
  edgeOnDASSize = DAS_size[edgeOn]


  faceOnDASMr = DAS_Mr[faceOn]
  midInclsDASMr = DAS_Mr[midincl]
  edgeOnInclsDASMr = DAS_Mr[edgeOn]
  
  
  graph11 = plot.Plots()
  FOsize_Data = plot.GraphData((faceOnMr, faceOnSize), 'b', 'step', False, 20, 'CALIFA b/a > 0.7')
  MIsize_Data = plot.GraphData((midInclsMr, midInclsSize), 'k', 'step', False, 20, 'CALIFA b/a between 0.3 and 0.7')
  EOsize_Data = plot.GraphData((edgeOnInclsMr, edgeOnSize), 'r', 'step', False, 20, 'CALIFA b/a < 0.3')

  DAS_FOsize_Data = plot.GraphData((faceOnDASMr, faceOnDASSize), DASColour, 'step', False, 20, 'DAS b/a > 0.7')
  DAS_MIsize_Data = plot.GraphData((midInclsDASMr, midInclsDASSize), DRASColour, 'step', False, 20, 'DAS b/a between 0.3 and 0.7')
  DAS_EOsize_Data = plot.GraphData((edgeOnInclsDASMr, edgeOnDASSize), 'grey', 'step', False, 20, 'DAS b/a < 0.3')


  bins = graph11.plotScatter([DAS_FOsize_Data, DAS_MIsize_Data, DAS_EOsize_Data, FOsize_Data, MIsize_Data, EOsize_Data], 'iso_size_vs_absmag', plot.PlotTitles("Physical size vs. absolute magnitude",  'M_r, mag', 'Isophotal size, kpc'), (-16, -24, 0, 90))

  exit()

  #physical size vs. absmag distribution
  size = convert(db.dbUtils.getFromDB('size', 'main', 'main'))
  absmag = convert(db.dbUtils.getFromDB('absmag', 'main', 'main'))
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'main', 'main'))/convert(db.dbUtils.getFromDB('isoA_r', 'main', 'main')) 

  
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnSize = size[faceOn]
  midInclsSize = size[midincl]
  edgeOnSize = size[edgeOn]


  faceOnMr = absmag[faceOn]
  midInclsMr = absmag[midincl]
  edgeOnInclsMr = absmag[edgeOn]
  
  
  
  graph11 = plot.Plots()
  FOsize_Data = plot.GraphData((faceOnMr, faceOnSize), 'b', 'step', False, 20, 'b/a > 0.7')
  MIsize_Data = plot.GraphData((midInclsMr, midInclsSize), 'k', 'step', False, 20, 'b/a between 0.3 and 0.7')
  EOsize_Data = plot.GraphData((edgeOnInclsMr, edgeOnSize), 'r', 'step', False, 20, 'b/a < 0.3')

  bins = graph11.plotScatter([FOsize_Data, MIsize_Data, EOsize_Data], 'reff_size_vs_absmag', plot.PlotTitles("Physical size vs. absolute magnitue",  'M_r, mag', 'R_e, kpc'), (-17.5, -24.5, 0, 35))

  exit()
  

  #plotting axis ratio distribution for different z with different M_r cuts
  isoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample'))
  graph11 = plot.Plots()
  isoA_rData = plot.GraphData((z, isoA_r), 'k', 'step', False, 20, 'Major isophotal axis')
  bins = graph11.plotScatter([isoA_rData], 'size_vs_z', plot.PlotTitles("Major axis vs. redshift",  'z', 'isoA_r'), (0.0045, 0.0305, 40, 85))
  exit()





  exit()
  #plotting axis ratio distribution for faint sample (M_r > 20.5)
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5')) 
  InclObs = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag > -20.5')) 
  graph11 = plot.Plots()

  #AxisRatioData = plot.GraphData((Incl), 'k', 'bar', True,'CALIFA mother sample, M_r > -20.5')
  AxisRatioData = plot.GraphData((Incl), 'k', 'bar', False, 20, 'CALIFA mother sample, M_r > -20.5')
  #ObsAxisRatioData = plot.GraphData((InclObs), 'r', 'bar', True,'Observed galaxies, M_r > -20.5')
  ObsAxisRatioData = plot.GraphData((InclObs), 'r', 'bar', False, 20, 'Observed galaxies, M_r > -20.5')
  bins = graph11.plotHist([AxisRatioData, ObsAxisRatioData], 'axis_ratio_faint', plot.PlotTitles("Axis ratio b/a distribution when M_r > -20.5",  'b/a', 'n'))


  #plotting axis ratio distribution for bright sample (M_r < 20.5)
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5')) 
  InclObs = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag < -20.5')) 
  graph11 = plot.Plots()

  #AxisRatioData = plot.GraphData((Incl), 'k', 'bar', True,'CALIFA mother sample, M_r < -20.5')
  #ObsAxisRatioData = plot.GraphData((InclObs), 'r', 'bar', True,'Observed galaxies, M_r < -20.5')
  AxisRatioData = plot.GraphData((Incl), 'k', 'bar', False, 20, 'CALIFA mother sample, M_r < -20.5')
  ObsAxisRatioData = plot.GraphData((InclObs), 'r', 'bar', False, 20, 'Observed galaxies, M_r < -20.5')
  
  bins = graph11.plotHist([AxisRatioData, ObsAxisRatioData], 'axis_ratio_bright', plot.PlotTitles("Axis ratio b/a distribution when M_r < -20.5",  'b/a', 'n'))

  exit()  
  #colour-sb-size
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'main', 'main'))
  sb = convert(db.dbUtils.getFromDB('sb', 'main', 'main'))
  isoA_r = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'main', 'main'))
  isoB_r = 0.396*convert(db.dbUtils.getFromDB('isoB_r', 'main', 'main'))
  ba = isoB_r/isoA_r
  
  print isoA_r
  
  faceOn = np.where(isoA_r > 60)
  midincl = np.where((isoA_r > 48) & (isoA_r < 60))
  edgeOn = np.where((isoA_r < 48))
  
  faceOnR = ba[faceOn]
  midInclsR = ba[midincl]
  edgeOnInclsR = ba[edgeOn]


  faceOnsb = sb[faceOn]
  midInclssb = sb[midincl]
  edgeOnInclssb = sb[edgeOn]
  
  graph = plot.Plots()
  FOSBData = plot.GraphData((faceOnR, faceOnsb), 'b', 'bar', False, 20, 'Face-on galaxies')  
  MidSBData = plot.GraphData((midInclsR, midInclssb), 'k', 'bar', False, 20, 'Face-on galaxies')  
  EOSBData = plot.GraphData((edgeOnInclsR, edgeOnInclssb), 'r', 'bar', False, 20, 'Face-on galaxies')  

  
  graph.plotScatter([FOSBData, MidSBData, EOSBData], 'SB_vs_ba', plot.PlotTitles("Surface brightness vs. b/a",  'Isophotal b/a', 'r band surface brightness, mag/arcsec^2'), (0, 1, 17, 23))  
  
  exit()
  #plotting axis ratio distribution for MS, DAS and RAS
  
  isoA_r_DAS = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS'))
  isoA_r_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS', ' where isoA_r > 20'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS', ' where isoA_r > 20'))
  print isoA_r_RAS.shape 
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  
  graph = plot.Plots()
  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio), CALIFAColour, 'step', True, 25,  'Mothersample')  
  DAS_AxisRatioData = plot.GraphData((isoA_r_DAS), DASColour, 'step', True, 25, 'DAS iso25_r')  
  RAS_AxisRatioData = plot.GraphData((isoA_r_RAS), RASColour, 'step', True, 25, 'RAS iso25_r, D_isoA > 20 pix (~8")')
  graph.plotHist([isoA_r_AxisRatioData, DAS_AxisRatioData, RAS_AxisRatioData], 'DAS_RAS_mothersample_axis_ratio_distribution_psf', plot.PlotTitles("b/a distributions for the mothersample, DAS and RAS", "b/a", "N"))  

  exit()

  r = convert(db.dbUtils.getFromDB('petroMag_r', 'DAS', 'DAS', ' where isoA_r > 0'))
  IsoA_r = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS', ' where isoA_r > 0')) 
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DAS', ' where isoA_r > 0'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS', ' where isoA_r > 0')) 
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnGR = IsoA_r[faceOn]
  midInclsGR = IsoA_r[midincl]
  edgeOnInclsGR = IsoA_r[edgeOn]
  graph = plot.Plots()
  FOData = plot.GraphData((r[faceOn], faceOnGR), 'b', 'bar', False, 20, 'b/a > 0.7') 
  MIData = plot.GraphData((r[midincl], midInclsGR), 'r', 'bar', False, 20, '0.3 < b/a < 0.7') 
  EOData = plot.GraphData((r[edgeOn], edgeOnInclsGR), 'g', 'bar', False, 20, 'b/a < 0.3') 
  graph.plotScatter([FOData, MIData, EOData], 'DAS_ISO_vs_r', plot.PlotTitles("IsoA_r vs. Petrosian r magnitude for IsoA_r-selected sample", 'r magnitude, mag',  'Isophotal b/a, arcsec'), (9, 17.8, 30, 90))
  

  r = convert(db.dbUtils.getFromDB('petroMag_r', 'RAS', 'RAS', ' where isoA_r > 0'))
  IsoA_r = 0.396*convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS', ' where isoA_r > 0')) 
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS', ' where isoA_r > 0'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS', ' where isoA_r > 0')) 
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnGR = IsoA_r[faceOn]
  midInclsGR = IsoA_r[midincl]
  edgeOnInclsGR = IsoA_r[edgeOn]
  graph = plot.Plots()
  FOData = plot.GraphData((r[faceOn], faceOnGR), 'b', 'bar', False, 20, 'b/a > 0.7') 
  MIData = plot.GraphData((r[midincl], midInclsGR), 'r', 'bar', False, 20, '0.3 < b/a < 0.7') 
  EOData = plot.GraphData((r[edgeOn], edgeOnInclsGR), 'g', 'bar', False, 20, 'b/a < 0.3') 
  graph.plotScatter([FOData, MIData, EOData], 'RAS_ISO_vs_r', plot.PlotTitles("IsoA_r vs. Petrosian r magnitude for z-selected sample", 'r magnitude, mag',  'Isophotal b/a, arcsec'), (11.5, 17.8, 0, 150))
  
  

 
  r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  IsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample')) 
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample')) 
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  faceOnGR = IsoA_r[faceOn]
  midInclsGR = IsoA_r[midincl]
  edgeOnInclsGR = IsoA_r[edgeOn]
  graph = plot.Plots()
  FOData = plot.GraphData((r[faceOn], faceOnGR), 'b', 'bar', False, 20, 'b/a > 0.7') 
  MIData = plot.GraphData((r[midincl], midInclsGR), 'r', 'bar', False, 20, '0.3 < b/a < 0.7') 
  EOData = plot.GraphData((r[edgeOn], edgeOnInclsGR), 'g', 'bar', False, 20, 'b/a < 0.3') 
  graph.plotScatter([FOData, MIData, EOData], 'ISO_vs_r', plot.PlotTitles("IsoA_r vs. Petrosian r magnitude", 'r magnitude, mag',  'Isophotal b/a, arcsec'), (10.5, 16, 30, 90))
  

  petro_r = convert(db.dbUtils.getFromDB('petroMag_r', 'CALIFA.sqlite', 'mothersample'))
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample')) 
  IsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample')) 
  #Incl_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS', ' where isoA_r > 20'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS', ' where isoA_r > 20'))
  #isoA_r_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS', ' where isoA_r > 20'))
  R90 = convert(db.dbUtils.getFromDB('r90', 'CALIFA.sqlite', 'nadine'))
  nadineIncl = convert(db.dbUtils.getFromDB('ba', 'CALIFA.sqlite', 'nadine'))
  
  graph = plot.Plots()
  #z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample'))
  
  Mean = np.mean(Incl)
  print Mean
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  nMean = np.mean(nadineIncl)
  #print round(nMean, 2)
  NfaceOn = np.where(nadineIncl > 0.7)
  Nmidincl = np.where((nadineIncl > 0.3) & (nadineIncl < 0.7))
  NedgeOn = np.where((nadineIncl < 0.3))
  
  faceOnGR = IsoA_r[faceOn]
  midInclsGR = IsoA_r[midincl]
  edgeOnInclsGR = IsoA_r[edgeOn]
  
  faceOnR90 = R90[NfaceOn]
  midInclsR90 = R90[Nmidincl]
  edgeOnInclsR90 = R90[NedgeOn]
  
  
  print 'sdss:', faceOnGR.shape, midInclsGR.shape, edgeOnInclsGR.shape
  print 'nadine:', faceOnR90.shape, midInclsR90.shape, edgeOnInclsR90.shape
  #FullGRData = plot.GraphData((IsoA_r), 'grey', 'step', False, 'total distribution')
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'bar', False, 50, 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', 'bar', False, 20, '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'bar', False, 20, 'b/a < 0.3')
  
  
  dataWidth =  max(faceOnR90) - min(faceOnR90)
  
  faceOnR90Data = plot.GraphData((faceOnR90), RASColour, 'bar', False, 15, 'b/a > 0.7')
  nmidInclsGRData = plot.GraphData((midInclsR90), 'k', 'bar', False, 60, '0.3 < b/a < 0.7')
  edgeOnR90Data = plot.GraphData((edgeOnInclsR90), DRASColour, 'bar', False, 15, 'b/a < 0.3')
  
  
  
  #zdata = plot.GraphData((z), 'grey', 'step', False, 'Total number')
  sizeMeasures = plot.GraphData((IsoA_r[0:937], R90), 'grey', 'step', False,20, 'size measures')
  baMeasures = plot.GraphData((Incl[0:937], nadineIncl), 'grey', 'step', False, 20, 'b/a measures')
  bins = graph.plotHist([nmidInclsGRData,edgeOnR90Data, faceOnR90Data], 'axisRatioHist_Nadines', plot.PlotTitles("R90 from growth curve photometry",  'R90, arcsec', 'n'), (10, 90, 0, 100))
  bins = graph.plotHist([midInclsGRData, edgeOnData, faceOnGRData], 'axisRatioHist', plot.PlotTitles("isoA_r",  'isoA_r, arcsec', 'n'), (40, 80, 0, 80))
  bins = graph.plotScatter([sizeMeasures], 'SizeMeasures', plot.PlotTitles("IsoA_r",  'IsoA_r, arcsec', 'R90, arcsec'), (40, 80, 10, 80))
  bins = graph.plotScatter([baMeasures], 'BAMeasures', plot.PlotTitles("b/a measures",  'growth curve photometry b/a', 'b/a from SDSS'))
  #faceOnGR = z[faceOn]
  #midInclsGR = z[midincl]
  #edgeOnInclsGR = z[edgeOn]
  
  
  


  exit()

  
  
  
    #plotting axis ratio distribution for bright sample (absmag_r < 20.5)
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag < -20.5'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample', ' where absmag < -20.5'))

  obsIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag < -20.5'))
  spiralIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample', ' where absmag < -20.5'))
  
  graph5 = plot.Plots()
  
  obsMagData = plot.GraphData((obsIsoA_r, observedSpiralsIncl), 'r', 'step', False, 'Observed galaxies with M_r < -20.5')

  spiralMagData = plot.GraphData((spiralIsoA_r, spiralsIncl), CALIFAColour,  'step', False, 'Mothersample, M_r < -20.5')
  graph5.plotScatter([spiralMagData, obsMagData], 'incl_vs_size_bright', plot.PlotTitles("Axis ratio vs. isophotal major axis for galaxies with M_r < -20.5", 'isoA_r, arcsec', 'Axis ratio'), (40, 80, 0, 1))
  
  #plotting axis ratio distribution for faint sample (absmag_r > 20.5)
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag > -20.5'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample', ' where absmag > -20.5'))

  graph5 = plot.Plots()
  obsIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed', ' where absmag > -20.5'))
  spiralIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample', ' where absmag > -20.5'))
  obsMagData = plot.GraphData((obsIsoA_r, observedSpiralsIncl), 'r', 'step', False, 'Observed galaxies')
  spiralMagData = plot.GraphData((spiralIsoA_r, spiralsIncl), CALIFAColour,  'step', False, 'Mothersample')
  graph5.plotScatter([spiralMagData, obsMagData], 'incl_vs_size_faint', plot.PlotTitles("Axis ratio vs. isophotal major axis for galaxies with M_r < -20.5", 'isoA_r, arcsec', 'Axis ratio'), (40, 80, 0, 1))
  

  

  
  
  #plotting Nadine's m_r vs. SDSS Petro_mag_r
  graph = plot.Plots()
  
  rMag = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  
  
  growth_curve_Mag = convert(db.dbUtils.getFromDB('r_mag', 'CALIFA.sqlite', 'nadine'))
  
  isoA_r_mag_Data = plot.GraphData((rMag), 'k', 'step', False, 'SDSS iso25_r')
  growth_curve_axis_ratioData = plot.GraphData((growth_curve_Mag), 'r', 'step', False, 'Growth curve photometry')
  
  graph.plotHist([isoA_r_mag_Data, growth_curve_axis_ratioData], 'growth_curve_and_iso25_app_mag_difference', plot.PlotTitles("Axis ratio vs. m_r from SDSS iso25_r, and growth curve", "m_r", "b/a"), 20)

  #plotting axis ratio distribution for different z with different M_r cuts
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5')) 
  z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5'))
  graph11 = plot.Plots()
  print z.shape
  
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  print 'faceOn', faceOn, 'midincl', midincl, 'edgeon', edgeOn
  
  faceOnGR = z[faceOn]
  midInclsGR = z[midincl]
  edgeOnInclsGR = z[edgeOn]
  
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'step', False, 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', 'step', False, '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'step', False,'b/a < 0.3')
  zdata = plot.GraphData((z), 'grey', 'step', False, 'Total number')
  bins = graph11.plotHist([zdata, faceOnGRData, midInclsGRData, edgeOnData], 'incl_vs_z_faint', plot.PlotTitles("Axis ratio vs. redshift, abs. mag > -20.5",  'z', 'n'), 10)

  
  #plotting axis ratio distribution for different z with r < 20.5
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5')) 
  z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5'))
  graph11 = plot.Plots()
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  faceOnGR = z[faceOn]
  midInclsGR = z[midincl]
  edgeOnInclsGR = z[edgeOn]
  
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'step', False, 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', 'step', False, '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'step', False,'b/a < 0.3')
  zdata = plot.GraphData((z), 'grey', 'step', False, 'Total number')
  bins = graph11.plotHist([zdata, faceOnGRData, midInclsGRData, edgeOnData], 'incl_vs_z_bright', plot.PlotTitles("Axis ratio vs. redshift, petroMag_r < 20.5",  'z', 'n'), 10)


  #plotting the distribution of major axis vs. m_r:
  graph9 = plot.Plots()
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  spiralIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  spiralMag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'mothersample'))
  
  
  faceOn = np.where(spiralsIncl > 0.7)
  midincl = np.where((spiralsIncl > 0.3) & (spiralsIncl < 0.7))
  edgeOn = np.where((spiralsIncl < 0.3))
  
  faceOnIncls = spiralIsoA_r[faceOn]
  midIncls = spiralIsoA_r[midincl]
  edgeOnIncls = spiralIsoA_r[edgeOn]
  
  faceOnMag = spiralMag[faceOn]
  midInclsMag = spiralMag[midincl]
  edgeOnInclsMag = spiralMag[edgeOn]
  
  faceOnInclsData = plot.GraphData((faceOnIncls, faceOnMag), 'b', 'step', False, ('b/a > 0.7'))
  midInclsData = plot.GraphData((midIncls, midInclsMag), 'k', 'step', False, ('0.3 < b/a < 0.7'))
  edgeOnInclsData = plot.GraphData((edgeOnIncls, edgeOnInclsMag), 'r', 'step', False, ('b/a < 0.3'))
  
  #print ks_2samp(faceOnIncls, edgeOnIncls)
  graph9.plotScatter([faceOnInclsData, midInclsData, edgeOnInclsData], 'isoA_vs_abs_mag_scatter', plot.PlotTitles("Isophotal major axis vs. abs. magnitude", 'isoA_R, arcsec', 'm_r'), (40, 80, -24, -16 ))

  
  
  #plotting the distribution of major axis vs. m_r:
  graph9 = plot.Plots()
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  spiralIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  spiralMag = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  
  
  faceOn = np.where(spiralsIncl > 0.7)
  midincl = np.where((spiralsIncl > 0.3) & (spiralsIncl < 0.7))
  edgeOn = np.where((spiralsIncl < 0.3))
  
  faceOnIncls = spiralIsoA_r[faceOn]
  midIncls = spiralIsoA_r[midincl]
  edgeOnIncls = spiralIsoA_r[edgeOn]
  
  faceOnMag = spiralMag[faceOn]
  midInclsMag = spiralMag[midincl]
  edgeOnInclsMag = spiralMag[edgeOn]
  
  faceOnInclsData = plot.GraphData((faceOnIncls, faceOnMag), 'b', ('b/a > 0.7'))
  midInclsData = plot.GraphData((midIncls, midInclsMag), 'k', ('0.3 < b/a < 0.7'))
  edgeOnInclsData = plot.GraphData((edgeOnIncls, edgeOnInclsMag), 'r', ('b/a < 0.3'))
  
  #print ks_2samp(faceOnIncls, edgeOnIncls)
  graph9.plotScatter([faceOnInclsData, midInclsData, edgeOnInclsData], 'isoA_vs_app_mag_scatter', plot.PlotTitles("Isophotal major axis vs. app. magnitude", 'isoA_R, arcsec', 'm_r'), (40, 80, 10.5, 16 ))
  
  #plotting axis ratio vs. m_r for MS, RAS, DAS, Nadine's MS data
  graph = plot.Plots()
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  rMag = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  
  growth_curve_axis_ratio = convert(db.dbUtils.getFromDB('ba', 'CALIFA.sqlite', 'nadine'))
  growth_curve_Mag = convert(db.dbUtils.getFromDB('r_mag', 'CALIFA.sqlite', 'nadine'))
  
  AxisRatio_RAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RAS'))
  rmag_RAS = convert(db.dbUtils.getFromDB('petroMag_r', 'RAS', 'RAS'))

  AxisRatio_DAS = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DAS'))
  rmag_DAS = convert(db.dbUtils.getFromDB('petroMag_r', 'DAS', 'DAS'))
  #print rmag_DAS.shape
  isoA_r_mag_Data = plot.GraphData((rMag, isoA_r_axis_ratio), 'k', 'step', False, 'SDSS iso25_r')
  growth_curve_axis_ratioData = plot.GraphData((growth_curve_Mag, growth_curve_axis_ratio), 'r', 'step', False, 'Growth curve photometry')
  RAS_Data = plot.GraphData((rmag_RAS, AxisRatio_RAS), RASColour, 'step', False, 'RAS')
  DAS_Data = plot.GraphData((rmag_DAS, AxisRatio_DAS), DASColour, 'step', False, 'DAS')
  graph.plotScatter([RAS_Data, DAS_Data, isoA_r_mag_Data, growth_curve_axis_ratioData], 'growth_curve_and_iso25_app_mag', plot.PlotTitles("Axis ratio vs. m_r from SDSS iso25_r, and growth curve", "m_r", "b/a"))

  #plotting the galaxies n vs. inclination distribution
  
  graph6 = plot.Plots()
  sersics = np.genfromtxt('sersics.txt')[:, 2]
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  
  diskIndices = np.where(sersics < 1)
  midIndices = np.where((1 <= sersics) & (sersics < 2.5))
  midHighIndices = np.where((2.5 <= sersics) & (sersics < 6))
  #highIndices = np.where((sersics >= 4) & (sersics < 6))
  
  diskInls = spiralsIncl[diskIndices]
  midInls =spiralsIncl[midIndices]
  midHighInls =spiralsIncl[midHighIndices]
  #highInls = spiralsIncl[highIndices]
  
  diskSersicsData = plot.GraphData((diskInls), 'k', 'bar', False,' n < 1')
  midSersicsData = plot.GraphData((midInls), 'g', 'stepfilled', False, '1 < n < 2.5')
  midHighSersicsData = plot.GraphData((midHighInls), 'b',  'step', False,'2.5 < n')
  #HighSersicsData = plot.GraphData((highInls), 'r', 'n > 4')
  #fakeInclsData = plot.GraphData((createFakeInclinations()), 'grey', 'random distribution of disks with q = 0.15')
  
  graph6.plotHist([midSersicsData, midHighSersicsData,  diskSersicsData], 'incl_vs_n_hist', plot.PlotTitles("Axis ratio vs. Sersic index", 'n', 'b/a'), 10)
  


  #plotting the observed spirals/mothersample spirals distribution, corrected for finite thickness
  
  obsSpiralsAxisRatio = np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'ObservedSpirals'))))/np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'ObservedSpirals'))))
  spiralAxisRatio =  np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))))/np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'Spirals'))))
  
  #print (spiralAxisRatio**2 - 0.13**2)
  SpiralsCorrectedIncl = np.sqrt((abs(spiralAxisRatio**2 - 0.15**2)/(1 - 0.15**2)))
  obsSpiralsCorrectedIncl = np.sqrt((abs(obsSpiralsAxisRatio**2 - 0.15**2)/(1 - 0.15**2)))
  #fakeInclsData = plot.GraphData((createFakeInclinations()), 'grey', 'random distribution of disks with q = 0.15')
  #print obsSpiralsCorrectedIncl, SpiralsCorrectedIncl
  graph2 = plot.Plots()
  
    
  obsSpiralsCorrectedInclData = plot.GraphData((obsSpiralsCorrectedIncl), 'r', 'bar', False, 'Observed sample')
  SpiralsCorrectedInclData = plot.GraphData((SpiralsCorrectedIncl), 'k', 'bar', False, 'Mother sample')
  graph2.plotHist([obsSpiralsCorrectedInclData, SpiralsCorrectedInclData], 'CorrectedInclDistribution', plot.PlotTitles("Inclination distribution, corrected for thickness", "Inclination", "Normalised number"), 20)    
  
  
  #plotting axis ratio distribution for both magnitude cuts
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  AxisRatio_DRAS = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DRAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DRAS'))
  AxisRatio_RDAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RDAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RDAS'))
  isoA_r_DRAS = convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DRAS'))*0.396
  isoA_r_RDAS = convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RDAS'))*0.396
  
  graph5 = plot.Plots()
  obsIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'observed'))
  spiralIsoA_r = convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  
  obsMagData = plot.GraphData((obsIsoA_r, observedSpiralsIncl), 'k', 'step', False, 'Observed galaxies')
  spiralMagData = plot.GraphData((spiralIsoA_r, spiralsIncl), CALIFAColour,  'step', False, 'Mothersample')
  DRASData = plot.GraphData((isoA_r_DRAS, AxisRatio_DRAS), DRASColour, 'step', False, 'DRAS galaxies')
  RDASData = plot.GraphData((isoA_r_RDAS, AxisRatio_RDAS), RDASColour, 'step', False, 'RDAS galaxies')
  graph5.plotScatter([ DRASData, RDASData,spiralMagData, obsMagData], 'incl_vs_size_all', plot.PlotTitles("Axis ratio vs. isophotal major axis for MS, DRAS and RDAS", 'isoA_r, arcsec', 'Axis ratio'), (40, 85, 0, 1))
  

  



  #plotting axis ratio distribution for different z with M_r > 20.5
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5')) 
  z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample', ' where absmag > -20.5'))
  graph11 = plot.Plots()
  print z.shape
  
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  print faceOn, midincl, edgeOn
  
  faceOnGR = z[faceOn]
  midInclsGR = z[midincl]
  edgeOnInclsGR = z[edgeOn]
  
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'step', False, 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', 'step', False, '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'step', False,'b/a < 0.3')
  zdata = plot.GraphData((z), 'grey', 'step', False, 'Total number')
  bins = graph11.plotHist([zdata, faceOnGRData, midInclsGRData, edgeOnData], 'incl_vs_z_faint', plot.PlotTitles("Axis ratio vs. redshift, petroMag_r > 20.5",  'z', 'n'), 10)

  #plotting axis ratio distribution for different z with r < 20.5
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5')) 
  z = convert(db.dbUtils.getFromDB('z', 'CALIFA.sqlite', 'mothersample', ' where absmag < -20.5'))
  graph11 = plot.Plots()
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  faceOnGR = z[faceOn]
  midInclsGR = z[midincl]
  edgeOnInclsGR = z[edgeOn]
  
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'step', False, 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', 'step', False, '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'step', False,'b/a < 0.3')
  zdata = plot.GraphData((z), 'grey', 'step', False, 'Total number')
  bins = graph11.plotHist([zdata, faceOnGRData, midInclsGRData, edgeOnData], 'incl_vs_z_bright', plot.PlotTitles("Axis ratio vs. redshift, petroMag_r < 20.5",  'z', 'n'), 10)



  
  #plotting axis ratio distribution for MS, RDAS and DRAS
  isoA_r_DRAS = convert(db.dbUtils.getFromDB('isoB_r', 'DAS', 'DRAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'DAS', 'DRAS'))
  isoA_r_RDAS = convert(db.dbUtils.getFromDB('isoB_r', 'RAS', 'RDAS'))/convert(db.dbUtils.getFromDB('isoA_r', 'RAS', 'RDAS'))
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  
  graph = plot.Plots()
  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio), CALIFAColour, 'step', False, 'Mothersample')  
  DRAS_AxisRatioData = plot.GraphData((isoA_r_DRAS), DRASColour, 'step', False,  'DRAS iso25_r')  
  RDAS_AxisRatioData = plot.GraphData((isoA_r_RDAS), RDASColour, 'step', False, 'RDAS iso25_r')
  graph.plotHist([isoA_r_AxisRatioData, DRAS_AxisRatioData, RDAS_AxisRatioData], 'DRAS_RDAS_mothersample_axis_ratio_distribution', plot.PlotTitles("b/a distributions for the mothersample, DRAS and RDAS", "b/a", "N"), 20)  

 
  #plotting the isoA_r/BA_90/growth curve ba spirals distribution
  graph = plot.Plots()
  isoA_r_axis_ratio = convert(db.dbUtils.getFromDB('isoB_r', 'CALIFA.sqlite', 'mothersample'))/convert(db.dbUtils.getFromDB('isoA_r', 'CALIFA.sqlite', 'mothersample'))
  growth_curve_axis_ratio = convert(db.dbUtils.getFromDB('ba', 'CALIFA.sqlite', 'nadine'))
  
  isoA_r_AxisRatioData = plot.GraphData((isoA_r_axis_ratio), 'k', 'step', False, 'SDSS iso25_r')
  growth_curve_axis_ratioData = plot.GraphData((growth_curve_axis_ratio), 'r', 'step', False, 'Growth curve photometry')
  BA90 = np.genfromtxt('BA90_filled.txt')[:, 2]
  print np.where(BA90 <0), BA90[np.where(BA90 <0)]
  BA90Data = plot.GraphData((BA90[np.where(BA90 >0)]), 'g', 'bar', False, 'BA_90 from NSAtlas')
	
  graph.plotHist([isoA_r_AxisRatioData, growth_curve_axis_ratioData, BA90Data], 'growth_curve_and_iso25_r_axis_ratio_distribution', plot.PlotTitles("b/a distributions from SDSS iso25_r, NSAtlas BA_90 and growth curve", "b/a", "N"), 20)
    




  #graph11.plotScatter([InlsData], 'incl_vs_z', plot.PlotTitles("Inclination vs. redshift",  'z', 'Inclination, deg'))



  
  sampleAxisRatio = convert(sami_db.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  gmAxisRatio = convert(sami_db.getFromDB('isoB_r', 'grandmother', 'grandmotherFL'))/convert(sami_db.getFromDB('isoA_r', 'grandmother', 'grandmotherFL'))
  print db.dbUtils.count('mothersample', 'observed') 
  
  observedAxisRatio = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'observed'))
  
  graph = plot.Plots()
  SampleAxisRatioData = plot.GraphData((sampleAxisRatio), 'k', 'Mother sample')
  gmAxisRatioData = plot.GraphData((gmAxisRatio), '0.5', 'Grandmother')
  observedAxisRatioData = plot.GraphData((observedAxisRatio), 'r', 'Observed')
  graph.plotHist([SampleAxisRatioData, gmAxisRatioData, observedAxisRatioData], 'axisRatioDistribution', plot.PlotTitles("Observed, mother and grandmother samples axis ratio distribution", "b/a", "N"), 10)

  
  #print ks_2samp(np.ndarray.flatten(sampleAxisRatio), np.ndarray.flatten(observedAxisRatio))
  #print ks_2samp(np.ndarray.flatten(gmAxisRatio), np.ndarray.flatten(observedAxisRatio))
  #np.savetxt('dump.txt', sampleData)
  #print sampleData.shape
  #sampleAMData = convert(sami_db.getFromDB('absmag', 'main', 'califa'))
  #gmAMData = convert(sami_db.getFromDB('absmag', 'grandmother', 'gm'))
  #plotScatter([gmAMData, gmSizeData, sampleAMData, sampleSizeData], ['black',  'red'],  'RED_BLACK_Size_vs_M', 'Size vs. Mr', (r'$M_{r}$'), (r'$M_{r}$'), ['Red -- CALIFA sample, black -- grandmother sample'], (-23, -14, 0, 30))
  
 
  #plotting the observed spirals/mothersample spirals distribution
  #observedSpiralsIncl = np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'ObservedSpirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'ObservedSpirals'))))
  #spiralsIncl = np.degrees(np.arccos(convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'Spirals'))))
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'ObservedSpirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'ObservedSpirals'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'Spirals'))
  
  graph = plot.Plots()
  observedSpiralsAxisRatioData = plot.GraphData((observedSpiralsIncl), 'r', 'Observed sample')
  SpiralsAxisRatioData = plot.GraphData((spiralsIncl), 'k', 'Mother sample')
  graph.plotHist([observedSpiralsAxisRatioData, SpiralsAxisRatioData], 'SpiralAxisRatioDistribution', plot.PlotTitles("Observed and mother sample spiral galaxies axis ratio distribution",  "b/a", "Normalised number"), 10)
  
  #print ks_2samp(np.ndarray.flatten(observedSpiralsAxisRatio), np.ndarray.flatten(spiralsAxisRatio))
  #obsSpiralsAxisRatio = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'ObservedSpirals'))/convert(db.dbUtils.getFromDB('isoA_r', 'mothersample', 'ObservedSpirals'))
  #spiralAxisRatio =  convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'Spirals'))
 
  
  #plotting the distribution of inclination vs. M:
  
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'observed'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  
  graph3 = plot.Plots()
  obsAbsMag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'observed'))
  spiralAbsMag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'mothersample'))
  obsAbsMagData = plot.GraphData((obsAbsMag, observedSpiralsIncl), 'r', 'Observed sample')
  spiralAbsMagData = plot.GraphData((spiralAbsMag, spiralsIncl), 'k', 'Mothersample')
  graph3.plotScatter([spiralAbsMagData, obsAbsMagData], 'Incl_vs_absMag', plot.PlotTitles("Axis ratio vs. absolute r magnitude", 'Mr, mag', 'b/a'))
  
  #plotting the distribution of inclination vs. m_r:
  
  observedSpiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'observed'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'observed'))
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))

  #print convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals', ' where isoA_r < 45'))
  
  graph4 = plot.Plots()
  obsMag = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'observed'))
  spiralMag = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  obsMagData = plot.GraphData((obsMag, observedSpiralsIncl), 'r', 'Observed galaxies')
  spiralMagData = plot.GraphData((spiralMag, spiralsIncl), 'k', 'Mothersample')
  graph4.plotScatter([spiralMagData, obsMagData], 'Incl_vs_appMag', plot.PlotTitles("Axis ratio vs. apparent r magnitude", 'm_r, mag', 'Axis ratio'))
  
  #plotting the distribution of inclination vs. colour:
  

  
  
  #inclination histogram for different magnitude bins
  
  
  
  
  
  #diskInls = spiralsIncl[diskIndices]
  #midInls =spiralsIncl[midIndices]
  #midHighInls =spiralsIncl[midHighIndices]
  #highInls = spiralsIncl[highIndices]

 
  #getting Sersics from NSAtlas
  #obsSpiralRa = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'ObservedSpirals'))
  #obsSpiralDec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'ObservedSpirals'))
  #spiralRa = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'Spirals'))
  #spiralDec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'Spirals'))
  #spiralID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals'))
   
  #sersicIndices = np.empty((len(spiralRa), 4))
  #sersics = readAtlas.find_in_atlas('SERSIC_N',  spiralRa, spiralDec, 3, spiralID)
  #print indicesData
  
  #np.savetxt('sersics.txt', sersics)
  #sersics = np.genfromtxt('sersics_nonrounded.txt', )
  
  #goodSersics = np.where(sersics[:, 2] != 0)
  #goodSersics = np.where(sersics[:, 2] != 0)
  #print sersics[goodSersics]
  #filledSersics = np.empty((len(spiralRa), 4))
  #filledSersics[goodSersics] = sersics[goodSersics]
  #filledSersics[zeroSersics] = readAtlas.find_in_atlas('SERSIC_N',  sersics[zeroSersics][:,0], sersics[zeroSersics][:,1], 2, sersics[zeroSersics][:,3])
  #filledSersics[:, 2] = filledSersics[:, 2]
  #print filledSersics[np.where(filledSersics[:, 2] > 9)]
  #filledSersics[np.where(filledSersics[:, 2] > 9)] = 10
  #np.savetxt('sersics.txt', filledSersics)
  
  #ra = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'mothersample'))
  #dec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'mothersample'))
  #ID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'mothersample'))
  #sersicIndices = np.empty((len(spiralRa), 4))
  #sersics = readAtlas.find_in_atlas('SERSIC_N',  ra, dec, 3, ID)
  #filledSersics = np.empty((len(ra), 4))
  #sersics = np.genfromtxt('all_sersics.txt')
  #print sersics[337, 2], sersics[656, 2]
  #sersics[337, 2] = 10
  #sersics[656, 2] = 10
  #np.savetxt('all_sersics.txt', sersics)
  #zeroSersics = np.where(sersics[:, 2] == 0)
  #goodSersics = np.where(sersics[:, 2] != 0)
  #filledSersics[goodSersics] = sersics[goodSersics]
  #filledSersics[zeroSersics] = readAtlas.find_in_atlas('SERSIC_N',  sersics[zeroSersics][:,0], sersics[zeroSersics][:,1], 2, sersics[zeroSersics][:,3])
  #np.savetxt('sersics_filled.txt', filledSersics)
  
  
  
  #plotting barredness vs. inclination

  graph7 = plot.Plots()
  barredness = np.genfromtxt("morph_nohead.csv", dtype = 'string', delimiter = ',')[:, 11]
  califaId = np.genfromtxt("morph_nohead.csv",  delimiter = ',')[:, 0]
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  
  barredIndices = np.where(barredness == ' B')
  unsureIndices = np.where(barredness == ' AB')
  noBarIndices = np.where(barredness == ' A')
  
  spiralID = convert((db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals')))
  spiralIndices = np.where(califaId == list(spiralID))
  
  spirals = Incl[spiralIndices[0]]

  barredSpiralsIndices =  findOverlap(spiralIndices[0], barredIndices[0])
  unsureSpiralsIndices = findOverlap(spiralIndices[0], unsureIndices[0])
  noBarSpiralsIndices = findOverlap(spiralIndices[0], noBarIndices[0])
  
  
  
  barInls = spirals[barredSpiralsIndices]
  unsureInls = spirals[unsureSpiralsIndices]
  noBarInls = spirals[noBarSpiralsIndices]
  
  print len(barInls), len(noBarInls)
  
  
  BarData = plot.GraphData((barInls), 'r', 'barred')
  UnsureBarData = plot.GraphData((unsureInls), 'b', 'undefined')
  noBarData = plot.GraphData((noBarInls), 'k', 'no bar')
  
  graph7.plotHist([BarData, UnsureBarData, noBarData], 'incl_vs_barredness_hist', plot.PlotTitles("Axis ratio vs. presence of bar", 'b/a', 'n'), 10, (0, 1, 0, 4))
  #plotting inclination vs. merger state
  Incl = db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample')/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  graph8 = plot.Plots()
  merger = np.genfromtxt("morph_nohead.csv", dtype = 'string', delimiter = ',')[:, 14]
  califaId = np.genfromtxt("morph_nohead.csv",  delimiter = ',')[:, 0]
  spiralID = convert((db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals')))
  spiralIndices = np.where(califaId == list(spiralID))
  mergerIndices = np.where(merger == ' M')
  isolatedIndices = np.where(merger == ' I')
  
  mergingSpiralsIndices =  findOverlap(spiralIndices[0], mergerIndices[0])
  isolatedSpiralsIndices = findOverlap(spiralIndices[0], isolatedIndices[0])
  spirals = Incl[spiralIndices[0]]
  
  mergerIncls = spirals[mergingSpiralsIndices]
  isolInls = spirals[isolatedSpiralsIndices]
  
  print len(mergerIncls), len(isolInls)
  
  mergData = plot.GraphData((mergerIncls), 'r', 'merging')
  isolInlsData = plot.GraphData((isolInls), 'k', 'isolated')
  
  graph8.plotHist([mergData, isolInlsData], 'incl_vs_merger_hist', plot.PlotTitles("Axis ratio vs. merger state",  'b/a', 'no.'), 10)
    

  #BA50 from NSAtlas
  #ra = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'mothersample'))
  #dec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'mothersample'))
  #ID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'mothersample'))
   
  #BA_indices = np.empty((len(ra), 4))
  #BA_indices = readAtlas.find_in_atlas('BA50',  ra, dec, 3, ID)
  #np.savetxt('BA_indices.txt', BA_indices)
  #print indicesData
  
  
  #ra = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'mothersample'))
  #dec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'mothersample'))
  #ID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'mothersample'))
  #sersicIndices = np.empty((len(spiralRa), 4))
  #sersics = readAtlas.find_in_atlas('SERSIC_N',  ra, dec, 3, ID)
  #filledBA = np.empty((len(ra), 4))
  #BA = np.genfromtxt('BA_indices.txt')
  #print BA[337, 2], BA[656, 2]
  #BA[337, 2] = 10
  #BA[656, 2] = 10
  #zeroBA = np.where(BA[:, 2] == 0)
  #goodBA = np.where(BA[:, 2] != 0)
  #filledBA[goodBA] = BA[goodBA]
  #filledBA[zeroBA] = readAtlas.find_in_atlas('BA50',  BA[zeroBA][:,0], BA[zeroBA][:,1], 2, BA[zeroBA][:,3])
  #print np.where(filledBA == 0)  
  #np.savetxt('BA_filled.txt', filledBA)
  
  

  BA50 = np.genfromtxt('BA_filled.txt')[:, 2]
  califaId = np.genfromtxt("BA_filled.txt")[:, 3]
  #print max(BA50) 
  spiralID = convert((db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals')))
  spiralIndices = np.where(califaId == list(spiralID))
  
  spiralBA50 = BA50[spiralIndices[0]]
  #print spiralBA50
  BA50Data = plot.GraphData((spiralBA50), 'r', 'spirals')
  #fakeInclsData = plot.GraphData((createFakeInclinations()), 'grey', 'random distribution of disks with q = 0.15')
  graph8 = plot.Plots()
  graph8.plotHist([BA50Data], 'BA50_hist', plot.PlotTitles("BA_50", 'b/a', 'no.'), 10, (0.1, 1, 0, 3))
 

  #BA90 from NSAtlas
  #ra = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'mothersample'))
  #dec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'mothersample'))
  #ID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'mothersample'))
   
  #BA_indices = np.empty((len(ra), 4))
  #BA_indices = readAtlas.find_in_atlas('BA90',  ra, dec, 3, ID)
  #np.savetxt('BA90_indices.txt', BA_indices)
  
  #sersicIndices = np.empty((len(spiralRa), 4))
  #sersics = readAtlas.find_in_atlas('SERSIC_N',  ra, dec, 3, ID)
  #filledBA = np.empty((len(ra), 4))
  #BA = np.genfromtxt('BA90_indices.txt')
  #print BA[337, 2], BA[656, 2]
  #BA[337, 2] = 10
  #BA[656, 2] = 10
  #zeroBA = np.where(BA[:, 2] == 0)
  #goodBA = np.where(BA[:, 2] != 0)
  #filledBA[goodBA] = BA[goodBA]
  #filledBA[zeroBA] = readAtlas.find_in_atlas('BA90',  BA[zeroBA][:,0], BA[zeroBA][:,1], 2, BA[zeroBA][:,3])
  #print np.where(filledBA == 0)  
  #np.savetxt('BA90_filled.txt', filledBA)
   

  BA90 = np.genfromtxt('BA90_filled.txt')[:, 2]
  #print max(BA90), min(BA90)
  califaId = np.genfromtxt("BA90_filled.txt")[:, 3]
   
  spiralID = convert((db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals')))
  spiralIndices = np.where(califaId == list(spiralID))
  
  spiralBA90 = BA90[spiralIndices[0]]

  BA90Data = plot.GraphData((spiralBA90), 'r', 'spirals')
  
  graph8 = plot.Plots()
  graph8.plotHist([BA90Data], 'BA90_hist', plot.PlotTitles("BA_90", 'b/a', 'n'), 10, (0.13, 1, 0, 2))


  #SERSIC_BA from NSAtlas
  #ra = convert(db.dbUtils.getFromDB('ra', 'mothersample', 'mothersample'))
  #dec = convert(db.dbUtils.getFromDB('dec', 'mothersample', 'mothersample'))
  #ID = convert(db.dbUtils.getFromDB('califa_id', 'mothersample', 'mothersample'))
   
  #BA_indices = np.empty((len(ra), 4))
  #BA_indices = readAtlas.find_in_atlas('SERSIC_BA',  ra, dec, 3, ID)
  #np.savetxt('SERSIC_BA90_indices.txt', BA_indices)
  
  #filledBA = np.empty((len(ra), 4))
 
  #BA = np.genfromtxt('SERSIC_BA90_indices.txt')
  #BA[337, 2] = 10
  #BA[656, 2] = 10
  #zeroBA = np.where(BA[:, 2] == 0)
  #goodBA = np.where(BA[:, 2] != 0)
  #filledBA[goodBA] = BA[goodBA]
  #filledBA[zeroBA] = readAtlas.find_in_atlas('SERSIC_BA',  BA[zeroBA][:,0], BA[zeroBA][:,1], 2, BA[zeroBA][:,3])
  #filledBA[337, 2] = 0
  #filledBA[656, 2] = 0 
  #np.savetxt('Sersic_BA_filled.txt', filledBA)
    
  
  
  BA = np.genfromtxt('Sersic_BA_filled.txt')[:, 2]
  
  califaId = np.genfromtxt("Sersic_BA_filled.txt")[:, 3]
   
  spiralID = convert((db.dbUtils.getFromDB('califa_id', 'mothersample', 'Spirals')))
  spiralIndices = np.where(califaId == list(spiralID))
  
  spiralBA = BA[spiralIndices[0]]
  
  goodSpiralBA = spiralBA[np.where(spiralBA > 0)]
  
  spiralBAData = plot.GraphData((goodSpiralBA), 'k', 'spiral galaxies b/a from Sersic fit')
  #fakeInclsData = plot.GraphData((createFakeInclinations()), 'grey', 'random distribution of disks with q = 0.15')
  graph9 = plot.Plots()
  graph9.plotHist([spiralBAData], 'Sersic_BA_hist', plot.PlotTitles("b/a from Sersic fit", 'b/a', 'n'), 10)
  
  #absmag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'mothersample'))
  #gr = (convert(db.dbUtils.getFromDB('petroMag_g', 'mothersample', 'mothersample')) - convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample')))
  
  #mgr = np.empty((len(absmag), 2))
  #mgr[:, 0] = absmag[:, 0]
  #mgr[:, 1] = gr[:, 0]
  
  #np.savetxt('gr_colour.txt', mgr)

  #plotting the inclination-coded distribution of colour vs. magnitude
  graph10 = plot.Plots()
  spiralsIncl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'mothersample'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'mothersample'))
  absmag = convert(db.dbUtils.getFromDB('absmag', 'mothersample', 'mothersample'))
  spiralMr = convert(db.dbUtils.getFromDB('petroMag_r', 'mothersample', 'mothersample'))
  g = convert(db.dbUtils.getFromDB('petroMag_g', 'mothersample', 'mothersample'))
  gr_colour =  g-spiralMr
  
  faceOn = np.where(spiralsIncl > 0.7)
  midincl = np.where((spiralsIncl > 0.3) & (spiralsIncl < 0.7))
  edgeOn = np.where((spiralsIncl < 0.3))
  
  faceOnMr = absmag[faceOn]
  #print faceOnMr.shape
  midIncMr = absmag[midincl]
  edgeOnMr = absmag[edgeOn]
  
  faceOnGR = gr_colour[faceOn]
  midInclsGR = gr_colour[midincl]
  edgeOnInclsGR = gr_colour[edgeOn]
  
  #print ks_2samp(midInclsGR, edgeOnInclsGR)
  
  faceOnGRData = plot.GraphData((faceOnGR, faceOnMr), 'b', ('b/a > 0.7'))
  midInclsGRData = plot.GraphData((midInclsGR, midIncMr), 'k', ('0.3 < b/a < 0.7'))
  edgeOnInclsGRData = plot.GraphData((edgeOnInclsGR, edgeOnMr), 'r', ('b/a < 0.3'))
  graph10.plotScatter([faceOnGRData, midInclsGRData, edgeOnInclsGRData], 'CMD_vs_incl', plot.PlotTitles("Magnitude vs. colour", 'g-r, mag', 'M_r, mag'))

  
  #createFakeInclinations()
  
 
  Incl = convert(db.dbUtils.getFromDB('isoB_r', 'mothersample', 'Spirals'))/convert(sami_db.getFromDB('isoA_r', 'mothersample', 'Spirals')) 
  z = convert(db.dbUtils.getFromDB('z', 'mothersample', 'Spirals'))
  graph11 = plot.Plots()
  
  faceOn = np.where(Incl > 0.7)
  midincl = np.where((Incl > 0.3) & (Incl < 0.7))
  edgeOn = np.where((Incl < 0.3))
  
  
  faceOnGR = z[faceOn]
  midInclsGR = z[midincl]
  edgeOnInclsGR = z[edgeOn]
  
  faceOnGRData = plot.GraphData((faceOnGR), 'b', 'b/a > 0.7')
  midInclsGRData = plot.GraphData((midInclsGR), 'k', '0.3 < b/a < 0.7')
  edgeOnData = plot.GraphData((edgeOnInclsGR), 'r', 'b/a < 0.3')
  zdata = plot.GraphData((z), 'grey', 'Total number')
  bins = graph11.plotHist([zdata, faceOnGRData, midInclsGRData, edgeOnData], 'incl_vs_z', plot.PlotTitles("Axis ratio vs. redshift",  'z', 'n'), 10)
  


  #graph11.plotScatter([InlsData], 'incl_vs_z', plot.PlotTitles("Inclination vs. redshift",  'z', 'Inclination, deg'))

  '''
 
if __name__ == '__main__':
    main()




