# -*- coding: utf-8 -*-
import pyfits
import numpy as np
import utils
#print head
import collections

#print atlasData.field('IAUNAME')

def find_mag(ra, dec):
    r_band = 4
    atlasFile = '/media/46F4A27FF4A2713B_/work/nsa_v0_1_2.fits'
    atlas = pyfits.open(atlasFile)
    atlasData = atlas[1].data
#    print atlas[1].columns.names

    head = atlas[1].header
    #print atlasData.field('RA')[15], atlasData.field('DEC')[15]
    for i in range(0, atlasData.shape[0]):
      if (round(atlasData.field('RA')[i], 3) == round(ra, 3)) and (round(atlasData.field('DEC')[i], 3) == round(dec, 3)):
	print 'i', i, round(atlasData.field('RA')[i], 3), round(atlasData.field('DEC')[i], 3)
	return utils.nmgy2mag(atlasData.field('NMGY')[i][r_band])
      #else:
	#print 'object not found!'

def find_in_atlas(arg, ra, dec):

    atlasFile = '../data/NSAtlas/nsa_v0_1_2.fits'
    atlas = pyfits.open(atlasFile)
    atlasData = atlas[1].data
#    print atlas[1].columns.names
    
    head = atlas[1].header
    #print atlasData.field('RA')[15], atlasData.field('DEC')[15]
    '''for i in range(0, atlasData.shape[0]):
      if (round(atlasData.field('RA')[i], 3) == round(ra, 3)) and (round(atlasData.field('DEC')[i], 3) == round(dec, 3)):
	print 'i', i, round(atlasData.field('RA')[i], 3), round(atlasData.field('DEC')[i], 3)
	return atlasData.field(arg)[i]
      #else:
	#print '''
    decs = np.where(np.round(atlasData.field('DEC'), 3) == round(dec, 3))
    ras = np.where(np.round(atlasData.field('RA'), 3) == round(ra, 3))
    print round(ra, 3), round(dec, 3)
    print ras[0]
    print decs[0]
    
    
    #if list(ras) == list(decs):
    
    a_multiset = collections.Counter(list(ras[0]))
    b_multiset = collections.Counter(list(decs[0]))
    overlap = list((a_multiset & b_multiset).elements())
    if overlap != []:
      ret =  atlasData.field(arg)[overlap]
    else:
      ret = 0
    print ret, 'returning '
    atlas.close()	
    return ret
    #print np.where((np.round(atlasData.field('RA'), 3) == round(ra, 3)))# and np.round(atlasData.field('DEC'), 3) == round(dec, 3)))
    #print np.where(np.round(atlasData.field('RA'), 3) == round(ra, 3))]#.field(arg)
