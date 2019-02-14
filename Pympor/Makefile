#########################################################################
## Create on compile program for Pympor utils                          ##
## Jun. 2015 by Guillaume THOMAS guillaume.thomas@astro.unistra.fr.    ##
## Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved. ##
#########################################################################


#--- gfortran -------------------------------
F90 = gfortran 
#--- gfortran -------------------------------



all: gravtoascii dumpdmparts
		chmod +x Pympor.py
				

gravtoascii:
		f2py -c --fcompiler=$(F90) -m gravtoascii ./data/gravtoascii.f90

dumpdmparts:
		f2py -c --fcompiler=$(F90) -m dumpdmparts ./data/dumpdmparts.f90

