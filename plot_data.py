import numpy as np
import db as db
import utils
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from scipy import integrate

#constants
c = 299792.458 
pi = 3.14159265
#     Cosmological parameters
H0 = 70.0 #km/s/Mpc, 1 Mpc= 3.08568025*1e19 km
Wm = 0.272
Wl = 0.728
Wh = c / H0 
Wk = 1 - Wm - Wl 
tH = (3.08568025*1e19/H0)/(3600 * 365 * 24*10e9) #Hubble time in Gyr
skyArea = 41252.9612494193 #degrees      
dH = c/H0 #in Mpc


#utilities
def E_z(z): #time derivative of log(a(t)), used for integration
    return 1/(np.sqrt((Wm*(1+z)**3 + Wk*(1+z)**2 + Wl)))

def comovingDistance(z):
    Dc = dH * 1000* integrate.quad(E_z, 0, z)[0] #in kpc
    return Dc

def angular2physical(reff, z): #return physical effective diameter of the galaxy in kpc
    return (np.radians(2*reff/3600) *(comovingDistance(z)) / (1 + z)**2)


params = {'backend': 'ps',
         'axes.labelsize': 10,
         'text.fontsize': 10,
         'legend.fontsize': 10,
         'xtick.labelsize': 8,
         'ytick.labelsize': 8,
         'text.usetex': True, 
         'font': 'serif',
         'font.size': 16,
         'ylabel.fontsize': 20}

from astroML.plotting import hist

dbDir = '../db/'


l_ids = utils.convert(db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'lucie', ' where a25 > 0'))
lucie_ids = ''
#id_length = 0
for i in l_ids:
    lucie_ids = lucie_ids+","+str(int(i))
#      id_length+=1
lucie_ids = '('+lucie_ids[1:]+')'

print l_ids.shape

redshift = utils.convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z', ' where califa_id in '+ lucie_ids))

print redshift.shape

ba = utils.convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine', ' where califa_id in '+ lucie_ids))
ba = np.reshape(ba, (ba.shape[0], ))

absmag = utils.convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect', ' where califa_id in '+ lucie_ids))
print absmag.shape, 'absmag'
iso_size = utils.convert(db.dbUtils.getFromDB('a25', dbDir+'CALIFA.sqlite', 'lucie', ' where califa_id in '+ lucie_ids))
print iso_size.shape, redshift.shape
size = np.empty((iso_size.shape))
for i, iso_size in enumerate(iso_size):
	size[i,0] = angular2physical(iso_size, redshift[i,0])
print size.shape
#np.savetxt('lucie_sizes.txt', size)
print np.min(ba), np.max(ba)
from astroML.plotting import hist
fig = plt.figure()
ax = fig.add_subplot(111)
#cm = cm.get_cmap('RdYlBu')
cm = cm.get_cmap('jet_r')
#ba = ba.tolist()
p = plt.scatter(absmag, size, c=ba,  marker='o', s=40,  cmap=cm, edgecolor="black", alpha=0.9) 
ax.axis([-17.8, -24, 0, 160])
minor_locator = plt.MultipleLocator(10)
Ymajor_locator = plt.MultipleLocator(50)  
major_locator = plt.MultipleLocator(1)      
Xminor_locator = plt.MultipleLocator(0.5)   
ax.xaxis.set_major_locator(major_locator)
ax.xaxis.set_minor_locator(Xminor_locator)     
ax.yaxis.set_major_locator(Ymajor_locator)
ax.yaxis.set_minor_locator(minor_locator)
#ax.set_ticks()
cb = fig.colorbar(p)
cb.set_label('b/a')
#cbar.set_ticks([0, 0.5, 1])
#cbar.ax.set_yticklabels(['0', '0.5', '1'])# vertically oriented colorbar
plt.xlabel("Absolute r magnitude, growth curve measurements")
plt.ylabel("Linear isophotal size from ellipse fits, kpc")

plt.savefig('sl_rel')
#SDSS sizes
size = np.genfromtxt('size.txt')
absmag = utils.convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect'))
fig = plt.figure()
ax = fig.add_subplot(111)
#cm = cm.get_cmap('RdYlBu')
#cm = cm.get_cmap('jet_r')
ba = utils.convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine'))
ba = np.reshape(ba, (ba.shape[0], ))

p = plt.scatter(absmag, size, c=ba,  marker='o', s=40,  cmap=cm, edgecolor="black", alpha=0.9)
ax.axis([-15, -26, 40, 80])

#minor_locator = plt.MultipleLocator(10)
#Ymajor_locator = plt.MultipleLocator(50)
#major_locator = plt.MultipleLocator(1)
#Xminor_locator = plt.MultipleLocator(0.5)
#ax.xaxis.set_major_locator(major_locator)
#ax.xaxis.set_minor_locator(Xminor_locator)
#ax.yaxis.set_major_locator(Ymajor_locator)
#ax.yaxis.set_minor_locator(minor_locator)
#ax.set_ticks()
cb = fig.colorbar(p)
cb.set_label('b/a')
#cbar.set_ticks([0, 0.5, 1])
#cbar.ax.set_yticklabels(['0', '0.5', '1'])# vertically oriented colorbar

plt.xlabel("Absolute r magnitude, growth curve measurements")
plt.ylabel("Linear isophotal size isoA_r, kpc")

plt.savefig('sl_rel_SDSS_size')


fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(ba, normed=False, color='grey', alpha=0.9)
plt.xlabel("b/a ratio")
plt.savefig('ba_hist_default')

full_redshift = utils.convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z'))
full_redshift = np.reshape(full_redshift, (full_redshift.shape[0], ))

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(full_redshift, bins='knuth', normed=False, color='grey', alpha=0.9)
plt.xlabel("redshift")
plt.savefig('z_hist')

#stellar masses
st_mass = utils.convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect'))
st_mass = np.reshape(st_mass, (st_mass.shape[0], ))
st_mass = np.log10(st_mass)
print st_mass
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(st_mass, bins='knuth', normed=False, color='grey', alpha=0.9)
plt.xlabel("log(Mst)")
plt.savefig('st_mass_hist')




































exit()
faint_ba = utils.convert(db.dbUtils.getFromDB('isoB_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where absmag < -19.5'))/utils.convert(db.dbUtils.getFromDB('isoA_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where absmag < -19.5'))
faint_ba = np.reshape(faint_ba, (faint_ba.shape[0], ))
hist(faint_ba, bins='knuth', normed=True)

ax = fig.add_subplot(122)
hist(ba, bins='blocks', normed=True, histtype='step')

hist(faint_ba, bins='blocks', normed=True, histtype='step')









'''
#SDSS
sdss_ba = utils.convert(db.dbUtils.getFromDB('isoB_r', dbDir+'RAS.sqlite', 'RAS'))/utils.convert(db.dbUtils.getFromDB('isoA_r', dbDir+'RAS.sqlite', 'RAS'))
sdss_ba = np.reshape(sdss_ba, (sdss_ba.shape[0], ))
fig = plt.figure(figsize=(20, 20))


ax = fig.add_subplot(121)
hist(sdss_ba, bins='knuth', normed=True, histtype='step')

ax = fig.add_subplot(122)
hist(sdss_ba, bins='blocks', normed=True, histtype='step')
'''


plt.show()
