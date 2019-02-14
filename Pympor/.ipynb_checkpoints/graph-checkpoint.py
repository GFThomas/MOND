#!/usr/bin/python
#########################################################################
##                Utils for plot the POR output                        ##
## Fev. 2015 by Guillaume THOMAS guillaume.thomas@astro.unistra.fr	  ##
## Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved. ##
#########################################################################

import numpy as np
from StringIO import StringIO
import math
import Pympor
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt
import random
import os


#######################
##      PLOTS        ##
####################### 

#----------------------
# Plot particules in 3D
#---------------------- 
def plot_part3D(no_output,data,time):
	title="%.2f Gyr" %(time/1000.0)
	fig=plt.figure()
	axes = fig.gca(projection='3d')
	axes.set_xlabel('x (kpc)')
	axes.set_ylabel('y (kpc)')
	axes.set_zlabel('z (kpc)')
	axes.scatter(data[:,0],data[:,1],data[:,2],s=0.1,alpha=0.5,color="blue")
	axes.set_xlim(-100,100)
	axes.set_ylim(-100,100)
	axes.set_zlim(-100,100)
	plt.title(title)
	plt.show("img_3D_%05d.png" %no_output)


#--------------------------------------
# Plot particules in the plan xy and xz 
#--------------------------------------
def plot_xy_xz(no_output,data,time):

	title="time = %.2f Gyr" %(time/1000.0)
	fig, (axes1, axes2) = plt.subplots(1, 2, figsize=(16,8))
	plt.title(title, fontweight='semibold', fontsize=14)

	axes1.set_xlabel('x (kpc)')
	axes1.set_ylabel('z (kpc)')
	axes1.annotate('Plan XZ', xy=(2, 1), xytext=(-95, 92))
	axes1.scatter(data[:,0],data[:,2],s=0.1,alpha=0.15,color="blue")	
	axes1.set_xlim(-100,100)
	axes1.set_ylim(-100,100)

	axes2.set_xlabel('x (kpc)')
	axes2.set_ylabel('y (kpc)')
	axes2.annotate('Plan XY', xy=(2, 1), xytext=(-95, 92))
	axes2.scatter(data[:,0],data[:,1],s=0.1,alpha=0.15,color="blue")
	axes2.set_xlim(-100,100)
	axes2.set_ylim(-100,100)

	fig.tight_layout()
	plt.show()
	#plt.savefig("img_%05d.png" %no_output)
	plt.close()



#######################
##       MAIN        ##
#######################
if __name__=='__main__':

	# Selection of the input files
	i_min=raw_input("No output min : ")
	i_min=int(i_min)
	i_max=raw_input("No output max : ")
	i_max=int(i_max)
	choix=raw_input("[0]:NEWTON  [1]:MOND ")
	choix=int(choix)

	for i in range (i_min,i_max+1):
		[data,time]=Pympor.chargement_part(i) # Particles loading 
		[part_MW,part_gal]=Pympor.filters(data,99000.0) # Selection of the kind of particles
		grav=Pympor.chargement_grav(i,choix) # Grav file loading
	
		plot_xy_xz(i,part_gal,time)
			
