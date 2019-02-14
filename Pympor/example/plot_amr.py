#!/usr/bin/python
#########################################################################
##                       Plot the AMR grid.                            ##
## May. 2015 by Guillaume THOMAS guillaume.thomas@astro.unistra.fr	  ##
## Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved. ##
#########################################################################

import numpy as np
import math
import matplotlib.pyplot as plt
import Pympor
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




def amr(data,plan,valeur_plan,valeur_min,valeur_max):
	
	# Selection of the plan
	if(plan==0): # Plan XY
		axe_a=data[:,0];axe_b=data[:,1];axe_c=data[:,2]
	elif(plan==1): # Plan XZ
		axe_a=data[:,0];axe_b=data[:,2];axe_c=data[:,1]
	elif(plan==2): # Plan YZ
		axe_a=data[:,1];axe_b=data[:,2];axe_c=data[:,0]
	
	# Selection of grid in the plan and in the area selection
	data=data[(((axe_c[:]+0.5*data[:,3])>=valeur_plan)*((axe_c[:]-0.5*data[:,3])<valeur_plan)*(axe_a[:]-0.5*data[:,3]>=valeur_min)*(axe_a[:]+0.5*data[:,3]<=valeur_max)*(axe_b[:]-0.5*data[:,3]>=valeur_min)*(axe_b[:]+0.5*data[:,3]<=valeur_max))]

	# Fill array with the good grid
	if(plan==0): # Plan XY
		axe_a=data[:,0];axe_b=data[:,1];axe_c=data[:,2]
	elif(plan==1): # Plan XZ
		axe_a=data[:,0];axe_b=data[:,2];axe_c=data[:,1]
	elif(plan==2): # Plan YZ
		axe_a=data[:,1];axe_b=data[:,2];axe_c=data[:,0]



	# Calcul of the matrix for AMR value
	inter_min=data[:,3].min()
	nb_interval=int((valeur_max-(valeur_min))/inter_min)
	indice=np.zeros((nb_interval+1,nb_interval+1))
	indmax=data[:,8].max()
	indmin=data[:,8].min()
	for ii in range (0,len(data[:,0])):
		diff=int(indmax-data[ii,8])
		for j in range (0,int(2**diff)):
			for k in range (0,int(2**diff)):
				a=axe_a[ii]-(2.0**diff-1.0)*(data[ii,3]/2.0**(diff+1.0))+j*(data[ii,3]/2.0**(diff))
				b=axe_b[ii]-(2.0**diff-1.0)*(data[ii,3]/2.0**(diff+1.0))+k*(data[ii,3]/2.0**(diff))	
				ind_a=int((a-valeur_min)/inter_min)
				ind_b=int((b-valeur_min)/inter_min)
				indice[ind_b,ind_a]=data[ii,8]

	##### PLOT AMR #####
	plt.imshow(indice, vmin=indmin, vmax=indmax, origin='lower',interpolation='nearest',
	           extent=[valeur_min, valeur_max, valeur_min, valeur_max])
	plt.xlabel('X (kpc)')
	plt.ylabel('Z (kpc)')	
	plt.title('AMR grid')
	cmap=plt.colorbar()
	cmap.set_label("Indice")
	labels = np.arange(indmin,indmax+1,1)
	loc    = labels
	cmap.set_ticks(loc)
	cmap.set_ticklabels(labels)
	plt.show()
	plt.close()
	####################






#######################
##       MAIN        ##
#######################
if __name__=='__main__':

	out=int(raw_input("No output (for the example enter 2): "))
	mond=int(raw_input("[0]:NEWTON or [1]:MOND (for the example enter 1): "))
	# 0 if it is a newtonian simulation (with Ramses) and 1 if MONDian simulation (with PoR)
	plan=int(raw_input("Which plan do you want to plot [0]:XY or [1]:XZ or [2]:YZ  (for the example enter 1): "))
	valeur_plan=float(raw_input("Which constant value for the other axe ? (for the example enter 7.0): "))
	print "Selection of a square area :"
	valeur_min=float(raw_input("Minimum value of the area (for the example enter -30): "))
	valeur_max=float(raw_input("Maximum value of the area (for the example enter 30): "))


	
	# Load the cells from the PoR output
	grid=Pympor.chargement_grav(out,mond) # Grid loading
# grid contain the information of the each cell of the AMR grid used by Ramses/PoR. Each row correspond to a cell with these different column  (x,y,z,inter,pot,ax,ay,az,ilevel,pot_N,ax_N,ay_N,az_N)
# These column correspond to :
# x: position of the center of cell along the x-axis (in kpc)
# y: position of the center of cell along the y-axis (in kpc)
# z: position of the center of cell along the z-axis(in kpc)
# inter: length of the cell (kpc)
# pot: potential at the center of the cell (in m2.s-2)
# ax: acceleration along the x-axis at the center of the cell (m.s-2)
# ay: acceleration along the y-axis at the center of the cell (m.s-2)
# az: acceleration along the z-axis at the center of the cell (m.s-2)
# ilevel: level of refinement of the grid
# pot_N: newtonian potential at the center of the cell (in m2.s-2)
# ax_N: newtonian acceleration along the x-axis at the center of the cell (m.s-2)
# ay_N: newtonian acceleration along the y-axis at the center of the cell (m.s-2)
# az_N: newtonian acceleration along the z-axis at the center of the cell (m.s-2)

## N.B. in case of a Newtonian output (mond=0), pot=pot_N, ax=ax_N, ...

	amr(grid,plan,valeur_plan,valeur_min,valeur_max) # Plot the amr grid

