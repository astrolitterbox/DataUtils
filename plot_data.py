import numpy as np
import db as db
from o_utils import convert as convert
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
         'axes.labelsize': 12,
         'text.fontsize': 12,
         'legend.fontsize': 10,
         'xtick.labelsize': 10,
         'ytick.labelsize': 10,
         'text.usetex': True, 
         'font': 'serif',
         'font.size': 12,
         'ylabel.fontsize': 12}

from astroML.plotting import hist
dbDir = '../../dev/db/'

#Stellar masses -- Rosa and SDSS kcorrect


r_ids = convert(db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'rosa'))
rosa_ids = ''
id_length = 0
for i in r_ids:
      rosa_ids = rosa_ids+","+str(int(i))
      id_length+=1
rosa_ids = '('+rosa_ids[1:]+')'
sdss_st_mass = convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect_sdss',  ' where califa_id in '+rosa_ids))
sdss_st_mass = np.log10(sdss_st_mass)
sdss_st_mass = np.reshape(sdss_st_mass, (sdss_st_mass.shape[0], ))

absmag = convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect', ' where califa_id in '+rosa_ids))
absmag = np.reshape(absmag, (absmag.shape[0], ))

rosa_st_mass = convert(db.dbUtils.getFromDB('starlight_mass', dbDir+'CALIFA.sqlite', 'rosa'))
rosa_st_mass = np.reshape(rosa_st_mass, (rosa_st_mass.shape[0], ))


fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
#print absmag.shape
cm_jet = cm.get_cmap('hsv')
p = plt.scatter(sdss_st_mass, rosa_st_mass, c=absmag, cmap=cm_jet, marker='o')
print ax.get_ylim()[0]
ax.axis([8.5, 11.5, 8.5, 11.5])
ax.plot([8.5, 11.5], [8.5, 11.5], 'black', lw=0.5)

cb = fig.colorbar(p)
cb.set_label('M_r')
plt.xlabel('Kcorrect stellar mass, SDSS photometry, log($M_\odot$)')
plt.ylabel('STARLIGHT stellar mass, log($M_\odot$)')
plt.savefig('st_mass_rosa_sdss')

exit()




st_mass_j = convert(db.dbUtils.getFromDB('mstar', dbDir+'CALIFA.sqlite', 'jakobs'))

petroR_ids = convert(db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'nadine', ' where r90 > 0'))

petro_ids = ''
id_length = 0
for i in petroR_ids:
      petro_ids = petro_ids+","+str(int(i))
      id_length+=1
petro_ids = '('+petro_ids[1:]+')'

#comparison between Petrosian R50 and GC circular hlr
gc_hlr = 0.396*convert(db.dbUtils.getFromDB('hlr', dbDir+'CALIFA.sqlite', 'r_tot_circ', ' where califa_id in '+petro_ids))
petroR_50 = convert(db.dbUtils.getFromDB('petroR50_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where califa_id in '+petro_ids))
#iso_size = convert(db.dbUtils.getFromDB('isoA_r', dbDir+'CALIFA.sqlite', 'mothersample'))


#iso size vs. z

concentration = convert(db.dbUtils.getFromDB('re/r90', dbDir+'CALIFA.sqlite', 'nadine', ' where califa_id in '+petro_ids))
iso_size = convert(db.dbUtils.getFromDB('isoA_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where califa_id in '+petro_ids))
iso_size = np.reshape(iso_size, (iso_size.shape[0], ))
size = np.empty((iso_size.shape))
redshift = convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z', ' where califa_id in '+petro_ids))

for i, iso_size in enumerate(iso_size):
	size[i,] = angular2physical(iso_size, redshift[i,])
print size.shape

cm_jet = cm.get_cmap('hsv')
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

p = plt.scatter(redshift, np.log10(size), c=concentration, cmap=cm_jet, marker='o', edgecolor="None")
	
#ax.axis([9.5, 15.5, 9.5, 15.5])
#ax.plot([9.5, 15.5], [9.5, 15.5], 'black', lw=0.5)
cb = fig.colorbar(p)
cb.set_label('concentration index')
plt.xlabel('z, corrected for infalls')

plt.ylabel('isophotal size, kpc')
plt.savefig('iso_size_z')

exit()


concentration = convert(db.dbUtils.getFromDB('re/r90', dbDir+'CALIFA.sqlite', 'nadine', ' where califa_id in '+petro_ids))
cm_jet = cm.get_cmap('hsv')
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

p = plt.scatter(petroR_50, gc_hlr, c=concentration, cmap = cm_jet, marker='o', edgecolor="None")
	
ax.axis([0, 38, 0, 38])
ax.plot([0, 38], [0, 38], 'black', lw=0.5)

cb = fig.colorbar(p)
cb.set_label('concentration index')
plt.xlabel('Petrosian R50, arcsec')

plt.ylabel('GC circular hlr, arcsec')
plt.savefig('petroR50_circ_hlr')
exit()

#comparison between elliptical photometry and circular annuli
gc_r = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'r_tot'))
circ_mag = convert(db.dbUtils.getFromDB('circ_mag', dbDir+'CALIFA.sqlite', 'r_tot_circ'))
ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine'))

cm_jet = cm.get_cmap('hsv')
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)

p = plt.scatter(circ_mag, gc_r, c=ba, cmap = cm_jet, marker='o', edgecolor="None")
	
ax.axis([9.5, 15.5, 9.5, 15.5])
ax.plot([9.5, 15.5], [9.5, 15.5], 'black', lw=0.5)

cb = fig.colorbar(p)
cb.set_label('b/a')
plt.xlabel('r mag in circular aperture')

plt.ylabel('r mag in elliptical aperture')
plt.savefig('circular_elliptical_mag_comparison')

#comparison between SDSS photometry and GC photometry stellar masses
gc_r = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'r_tot'))
petroMag_r = convert(db.dbUtils.getFromDB('petroMag_r', dbDir+'CALIFA.sqlite', 'mothersample'))

st_mass = convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect'))
st_mass = np.log10(st_mass)
st_mass = np.reshape(st_mass, (st_mass.shape[0], ))
sdss_st_mass = convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect_sdss'))
sdss_st_mass = np.log10(sdss_st_mass)
sdss_st_mass = np.reshape(sdss_st_mass, (sdss_st_mass.shape[0], ))

cm_jet = cm.get_cmap('hsv')
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
mag_diff = petroMag_r-gc_r
mag_diff[np.where(mag_diff > 2)] = 0
p = plt.scatter(sdss_st_mass, st_mass, c=mag_diff, cmap = cm_jet, marker='o', edgecolor="None")
	
ax.axis([8.5, 11.5, 8.5, 11.5])
ax.plot([8.5, 11.5], [8.5, 11.5], 'black', lw=0.5)

cb = fig.colorbar(p)
cb.set_label('PetroMag_r - gc r magnitude')
plt.xlabel('SDSS photometry stellar mass, log($M_\odot$)')
plt.ylabel('GC photometry stellar mass, log($M_\odot$)')
plt.savefig('st_mass_SDSS_GC_comparison')



cm = cm.get_cmap('jet_r')

#a/b, axis ratio, compared with RAS
#ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine'))
ba = convert(db.dbUtils.getFromDB('isoB_r/isoA_r', dbDir+'CALIFA.sqlite', 'mothersample'))
ba = np.reshape(ba, (ba.shape[0], ))

ba_ras = convert(db.dbUtils.getFromDB('isoB_r/isoA_r', dbDir+'RAS.sqlite', 'RAS', ' where isoA_r > 15 and petroMag_r < 20'))
ba_ras = np.reshape(ba_ras, (ba_ras.shape[0], ))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
n, bins, patches = hist(ba, normed=True, bins='knuth', color='red', alpha=0)
hist(ba_ras, normed=True, bins=bins, label='RAS b/a', histtype='stepfilled', color='grey', alpha=1, hatch='o')
hist(ba, normed=True, bins=bins, alpha=0.8, label='CALIFA b/a', color='red')
plt.legend()
plt.xlabel("Isophotal b/a")
plt.savefig('ba_hist_RAS', bbox_inches='tight')
exit()
#apparent r magnitude comparison
concentration = convert(db.dbUtils.getFromDB('re/r90', dbDir+'CALIFA.sqlite', 'nadine', ' where califa_id in '+petro_ids))

r_mag = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'r_tot', ' where califa_id in '+petro_ids))
petroMag_r = convert(db.dbUtils.getFromDB('petromag_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where califa_id in'+petro_ids))
r_mag = np.reshape(r_mag, (r_mag.shape[0], ))
petroMag_r = np.reshape(petroMag_r, (petroMag_r.shape[0], ))

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
#hist(petroMag_r, bins='knuth', normed=False, color='grey', alpha=1, label='SDSS petroMag_r')
#hist(r_mag, bins='knuth', normed=False, color='red', alpha=0.5, label='GC r mag')
p = plt.scatter(petroMag_r, r_mag, s=10, c=concentration, cmap=cm, label='GC r mag vs. SDSS mag', edgecolor='None')
ax.axis([9.5, 16, 9.5, 16])
ax.plot([9.5, 16], [9.5, 16])
plt.legend()
plt.xlabel("petroMag_r, mag")
plt.ylabel("GC r mag, mag")
cb = fig.colorbar(p)
cb.set_label('r_50/r_90')
plt.savefig('apparent_mag_comparison_conc')

#apparent u magnitude


u_mag = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'u_tot', ' where califa_id in '+petro_ids))
petroMag_u = convert(db.dbUtils.getFromDB('petromag_u', dbDir+'CALIFA.sqlite', 'mothersample', ' where califa_id in'+petro_ids))
u_mag = np.reshape(u_mag, (u_mag.shape[0], ))
petroMag_u = np.reshape(petroMag_u, (petroMag_u.shape[0], ))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
p = plt.scatter(petroMag_u, u_mag, s=10, c=concentration, cmap=cm, label='GC u mag vs. SDSS mag', edgecolor='None')
ax.axis([12.5, 18, 12.5, 18])
ax.plot([12.5, 18], [12.5, 18], c='black')
plt.xlabel("petroMag_u, mag")
plt.ylabel("GC u mag, mag")
plt.legend()
cb = fig.colorbar(p)
cb.set_label('r_50/r_90')
plt.savefig('apparent_u_mag_comparison_conc')

#comparison with NSAtlas

r_mag = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'r_tot', ' where califa_id in '+petro_ids))
atlas_mag = convert(db.dbUtils.getFromDB('r_mag', dbDir+'CALIFA.sqlite', 'atlas', ' where califa_id in'+petro_ids))
r_mag = np.reshape(r_mag, (r_mag.shape[0], ))
atlas_mag = np.reshape(atlas_mag, (atlas_mag.shape[0], ))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
p = plt.scatter(atlas_mag, r_mag, s=10, c=concentration, cmap=cm, label='GC r mag vs. NSAtlas', edgecolor='None')
ax.axis([9.5, 15.5, 9.5, 15.5])
ax.plot([9.5, 15.5], [9.5, 15.5], c='black')
plt.xlabel("NASA SLOAN Atlas mag")
plt.ylabel("GC r mag, mag")
plt.legend()
cb = fig.colorbar(p)
cb.set_label('r_50/r_90')
plt.savefig('nsAtlas_comparison')



exit()
#CMD
g_mag = convert(db.dbUtils.getFromDB('el_mag', dbDir+'CALIFA.sqlite', 'g_tot'))
g_mag = np.reshape(g_mag, (g_mag.shape[0], ))
absmag = convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect'))
absmag = np.reshape(absmag, (absmag.shape[0], ))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
p = plt.scatter(absmag, r_mag - g_mag)
plt.savefig('cmd')
print g_mag.shape, absmag.shape




#absolute r mag

#absmag_r = convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect'))
#petroMag_r = convert(db.dbUtils.getFromDB('petromag_r', dbDir+'CALIFA.sqlite', 'mothersample'))

#absmag_r = -1*np.reshape(absmag_r, (absmag_r.shape[0], ))
#petroMag_r = np.reshape(petroMag_r, (petroMag_r.shape[0], ))
print absmag.shape
fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
#hist(petroMag_r, bins='knuth', normed=False, color='grey', alpha=1)
hist(absmag, bins=10, normed=False, color='red', alpha=0.8)

plt.xlabel("Absolute r magnitude, mag")
plt.savefig('absolute_mag_hist')



#Stellar mass

st_mass = convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect'))
#petroMag_r = convert(db.dbUtils.getFromDB('petromag_r', dbDir+'CALIFA.sqlite', 'mothersample'))
st_mass = np.log10(st_mass)
st_mass = np.reshape(st_mass, (st_mass.shape[0], ))
st_mass_j = np.log10(st_mass_j)
st_mass_j = np.reshape(st_mass_j, (st_mass_j.shape[0], ))
#petroMag_r = np.reshape(petroMag_r, (petroMag_r.shape[0], ))

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
#hist(petroMag_r, bins='knuth', normed=False, color='grey', alpha=1)
hist(st_mass, bins='knuth', normed=False, color='red', alpha=0.8)
hist(st_mass_j, bins='knuth', normed=False, color='grey', alpha=0.6)
plt.xlabel("Stellar mass, log($M_\odot$)")
plt.savefig('st_mass_hist')




#redshift distribution
redshift = convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z'))
redshift = np.reshape(redshift, (redshift.shape[0], ))

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
n, bins, patches = hist(redshift, bins='knuth', normed=False, color='grey', alpha=0.9)
np.savetxt('bins.txt', bins)
np.savetxt('n.txt', n)
plt.xlabel("Redshift")
plt.savefig('z_hist')

#a/b, axis ratio
ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine'))
ba = np.reshape(ba, (ba.shape[0], ))

fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(ba, normed=False, bins='knuth', color='grey', alpha=0.9)
plt.xlabel("b/a ratio")
plt.savefig('ba_hist')







r_ids = convert(db.dbUtils.getFromDB('califa_id', dbDir+'CALIFA.sqlite', 'rosa'))
rosa_ids = ''
id_length = 0
for i in r_ids:
      rosa_ids = rosa_ids+","+str(int(i))
      id_length+=1
rosa_ids = '('+rosa_ids[1:]+')'

#Stellar mass

st_mass = convert(db.dbUtils.getFromDB('st_mass', dbDir+'CALIFA.sqlite', 'kcorrect', ' where califa_id in '+rosa_ids))
#petroMag_r = convert(db.dbUtils.getFromDB('petromag_r', dbDir+'CALIFA.sqlite', 'mothersample'))
st_mass = np.log10(st_mass)
st_mass = np.reshape(st_mass, (st_mass.shape[0], ))
rosa_st_mass = convert(db.dbUtils.getFromDB('starlight_mass', dbDir+'CALIFA.sqlite', 'rosa'))
rosa_st_mass = np.reshape(rosa_st_mass, (rosa_st_mass.shape[0], ))

#st_mass_j = np.log10(st_mass_j)
#st_mass_j = np.reshape(st_mass_j, (st_mass_j.shape[0], ))
#petroMag_r = np.reshape(petroMag_r, (petroMag_r.shape[0], ))

fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
#hist(petroMag_r, bins='knuth', normed=False, color='grey', alpha=1)
n, bins, patches = hist(st_mass, bins='knuth', normed=False, color='red', alpha=0.8, label='kcorrect $M_{st}$')

hist(rosa_st_mass, bins=bins, normed=False, color='grey', alpha=0.6, label='STARLIGHT $M_{st}$')
plt.legend()
plt.xlabel("Stellar mass, log($M_\odot$)")
plt.savefig('st_mass_hist_rosa')

#stellar mass scatterplot
#cm = cm.get_cmap('jet')
iso_size = convert(db.dbUtils.getFromDB('isoA_r',  dbDir+'CALIFA.sqlite', 'mothersample', ' where califa_id in '+ rosa_ids))

absmag_rosa = convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect', ' where califa_id in '+ rosa_ids))
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
print absmag.shape
p = plt.scatter(st_mass, rosa_st_mass, c=absmag_rosa, cmap = cm, marker='o')
print ax.get_ylim()[0]
ax.axis([8.5, 11.5, 8.5, 11.5])
ax.plot([8.5, 11.5], [8.5, 11.5], 'black', lw=0.5)

cb = fig.colorbar(p)
cb.set_label('M_r')
plt.xlabel('Kcorrect stellar mass, log($M_\odot$)')
plt.ylabel('STARLIGHT stellar mass, log($M_\odot$)')
plt.savefig('st_mass_comparison')

#print l_ids.shape

#redshift = convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z', ' where califa_id in '+ lucie_ids))

#print redshift.shape

#ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine', ' where califa_id in '+ lucie_ids))
#ba = np.reshape(ba, (ba.shape[0], ))

#absmag = convert(db.dbUtils.getFromDB('r', dbDir+'CALIFA.sqlite', 'kcorrect', ' where califa_id in '+ lucie_ids))
#print absmag.shape, 'absmag'
#iso_size = convert(db.dbUtils.getFromDB('a25', dbDir+'CALIFA.sqlite', 'lucie', ' where califa_id in '+ lucie_ids))
#print iso_size.shape, redshift.shape
iso_size = convert(db.dbUtils.getFromDB('isoA_r', dbDir+'CALIFA.sqlite', 'mothersample'))
iso_size = np.reshape(iso_size, (iso_size.shape[0], ))
size = np.empty((iso_size.shape))
for i, iso_size in enumerate(iso_size):
	size[i,] = angular2physical(iso_size, redshift[i,])
print size.shape
#np.savetxt('lucie_sizes.txt', size)
print np.min(ba), np.max(ba)
from astroML.plotting import hist
fig = plt.figure()
ax = fig.add_subplot(111)
#cm = cm.get_cmap('RdYlBu')
#ba = ba.tolist()
p = plt.scatter(absmag, size, c=ba,  marker='o', s=40,  cmap=cm, edgecolor="black", alpha=0.9) 
ax.axis([-17, -24.5, 0, 90])
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
plt.ylabel("Linear isophotal size, kpc")

plt.savefig('sl_rel')
#ba = convert(db.dbUtils.getFromDB('ba', dbDir+'CALIFA.sqlite', 'nadine'))
#ba = np.reshape(ba, (ba.shape[0], ))
'''
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
'''

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(ba, normed=False, color='grey', alpha=0.9)
plt.xlabel("b/a ratio")
plt.savefig('ba_hist_default')

full_redshift = convert(db.dbUtils.getFromDB('z', dbDir+'CALIFA.sqlite', 'ned_z'))
full_redshift = np.reshape(full_redshift, (full_redshift.shape[0], ))

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15)
hist(full_redshift, bins='knuth', normed=False, color='grey', alpha=0.9)
plt.xlabel("redshift")
plt.savefig('z_hist')

































exit()
faint_ba = utils.convert(db.dbUtils.getFromDB('isoB_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where absmag < -19.5'))/utils.convert(db.dbUtils.getFromDB('isoA_r', dbDir+'CALIFA.sqlite', 'mothersample', ' where absmag < -19.5'))
faint_ba = np.reshape(faint_ba, (faint_ba.shape[0], ))
hist(faint_ba, bins='knuth', normed=True)

ax = fig.add_subplot(122)
hist(ba, bins='blocks', normed=True, histtype='step')

hist(faint_b/a, bins='blocks', normed=True, histtype='step')









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
