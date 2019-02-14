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
				ind_a=math.floor((a-valeur_min)/inter_min)
				ind_b=math.floor((b-valeur_min)/inter_min)
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

	out=raw_input("No output : ")
	out=int(out)
	mond=raw_input("[0]:NEWTON or [1]:MOND  ")
	mond=int(mond)
	plan=raw_input("[0]:XY or [1]:XZ or [2]:YZ  ")
	plan=int(plan)
	valeur_plan=raw_input("Which constant value for the other axe ? ")
	valeur_plan=float(valeur_plan)
	print "Selection of a square area :"
	valeur_min=raw_input("Minimum value of the area ")
	valeur_min=float(valeur_min)
	valeur_max=raw_input("Maximum value of the area ")
	valeur_max=float(valeur_max)


	data=Pympor.chargement_grav(out,mond) # Grid loading
	amr(data,plan,valeur_plan,valeur_min,valeur_max) # Plot the amr grid

