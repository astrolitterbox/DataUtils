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

def find_in_atlas(arg, ra, dec, roundParam, spiralID):

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
    ret = np.empty((len(ra), 4))
    for i in range(0, len(ra)):
      
      decs = np.where(np.round(atlasData.field('DEC'), roundParam) == round(dec[i], roundParam))
      ras = np.where(np.round(atlasData.field('RA'), roundParam) == round(ra[i], roundParam))
      #print round(ras, 3), round(decs, 3)
      print 'ra: ', ras[0], 'dec', decs[0]
    
      a_multiset = collections.Counter(list(ras[0]))
      b_multiset = collections.Counter(list(decs[0]))
      overlap = list((a_multiset & b_multiset).elements())
      print 'o', overlap, spiralID[i]
      if overlap != []:
	ret[i][0] =  atlasData.field('ra')[overlap]
	ret[i][1] = atlasData.field('dec')[overlap]      
	ret[i][2] =  atlasData.field(arg)[overlap]
	
      else:
	ret[i][0] =  ra[i]
	ret[i][1] = dec[i]
	ret[i][2] = 0
      ret[i][3] = spiralID[i]
      print 'sersic index ', ret[i]
    
    atlas.close()	
    return ret
    #print np.where((np.round(atlasData.field('RA'), 3) == round(ra, 3)))# and np.round(atlasData.field('DEC'), 3) == round(dec, 3)))
    #print np.where(np.round(atlasData.field('RA'), 3) == round(ra, 3))]#.field(arg)
