#!/usr/bin/python
#########################################################################
## Utils for read POR outputs                                          ##
## Fev. 2015 by Guillaume THOMAS guillaume.thomas@astro.unistra.fr	  ##
## Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved. ##
#########################################################################

import numpy as np
from StringIO import StringIO
import math
import dumpdmparts
import gravtoascii



####################################
######  Reading subroutines  #######
####################################
# Read the particles data of POR
def chargement_part(nb_output):
	name="output_%05d" %(nb_output)
	print "File :",name
	data = np.zeros((10000000,9), dtype='f', order='Fortran')
	data[:,0:8],time=dumpdmparts.star2list(name) #data(x,y,z,vx,vy,vz,mass,id) and time is the age of the simulation in Myr"
	data=data[:,0:7]
	print "part file read"
	return data,time


# Filter between the MW particles and the particles of the dwarf galaxy
def filters(data,m_threshold): # m_threshold in Solar mass
	part_MW=data[(data[:,6]>m_threshold)]
	part_gal= data[(data[:,6]<m_threshold)*(data[:,6]>0.0)]	
	return part_MW,part_gal


def chargement_grav(nb_output,mond): # mond = 1 if MONDian case, 0 if Newtonian case
	name="output_%05d" %(nb_output)
	data = np.zeros((100000000,13), dtype='f', order='Fortran')
	data=gravtoascii.grav2ascii(name,mond) #data(x,y,z,inter,pot,ax,ay,az,ilevel,pot_N,ax_N,ay_N,az_N) 
	print "grav file read"
	data_filtre=data[(data[:,8]>0.0)]
	return data_filtre
#####################################################################################




