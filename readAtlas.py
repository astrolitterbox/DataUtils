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

def find_in_atlas(arg, ra, dec, roundParam):

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
    ret = np.empty((len(ra), 3))
    for i in range(0, len(ra)):
      
      decs = np.where(np.round(atlasData.field('DEC'), roundParam) == round(dec[i], roundParam))
      ras = np.where(np.round(atlasData.field('RA'), roundParam) == round(ra[i], roundParam))
      #print round(ras, 3), round(decs, 3)
      print 'ra: ', ras[0], 'dec', decs[0], 'i', i
    
      a_multiset = collections.Counter(list(ras[0]))
      b_multiset = collections.Counter(list(decs[0]))
      overlap = list((a_multiset & b_multiset).elements())
      print 'o', overlap, len(overlap), roundParam
      if len(overlap) == 1:
	ret[i][0] =  atlasData.field('ra')[overlap]
	ret[i][1] = atlasData.field('dec')[overlap]      
	ret[i][2] =  round(atlasData.field(arg)[overlap], 1)
      elif overlap == []:
	try:
	  ret[i][0] =  ra[i]
	  ret[i][1] = dec[i]
	  print 'repeat', find_in_atlas(arg, ra[i], dec[i], roundParam-1)[:, 2]
	  ret[i][2] = find_in_atlas(arg, ra[i], dec[i], roundParam-1)[:, 2]
	except TypeError:
	  print 'not detected'
	  ret[i][0] =  ra[i]
	  ret[i][1] = dec[i]	
	  ret[i][2] = 0	  
      else:
	try:
	  ret[i][0] =  ra[i]
	  ret[i][1] = dec[i]
	  print 'repeat', find_in_atlas(arg, ra[i], dec[i], roundParam+1)[:, 2]
	  ret[i][2] = find_in_atlas(arg, ra[i], dec[i], roundParam+1)[:, 2]
	except TypeError:
	  ret[i][0] =  ra[i]
	  ret[i][1] = dec[i]
	  print 'too many values, ambiguous detection'
	  ret[i][2] = 0
      print 'sersic index ', ret[i]
    
    atlas.close()	
    return ret
    #print np.where((np.round(atlasData.field('RA'), 3) == round(ra, 3)))# and np.round(atlasData.field('DEC'), 3) == round(dec, 3)))
    #print np.where(np.round(atlasData.field('RA'), 3) == round(ra, 3))]#.field(arg)
