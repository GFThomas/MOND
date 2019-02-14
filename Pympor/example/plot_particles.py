#!/usr/bin/python
#########################################################################
##                Utils for plot the POR output                        ##
## Fev. 2015 by Guillaume THOMAS guillaume.thomas@astro.unistra.fr	  ##
## Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved. ##
#########################################################################


# Package needed to run the code
import numpy as np
from StringIO import StringIO
import math
import Pympor
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.pyplot as plt




n_output=int(raw_input("No output to read (for the example enter 2): ")) # Number of the output to read

[particles,time]=Pympor.chargement_part(n_output) # Particles loading

# particles contain the particles of the simulation with as column : x(kpc) y(kpc) z(kpc) vx(kpc) vy(kpc) vz(kpc) mass(M_sun)
# time correspond to the time of the output (in Myr)

print "Particles readed"
print "Number of particles :",len(particles)


# Plot the position of the particles
title="time = %.2f Gyr" %(time/1000.0)
fig, (axes1, axes2) = plt.subplots(1, 2, figsize=(10.5,5))
plt.title(title, fontweight='semibold', fontsize=14)

axes1.set_xlabel('x (kpc)')
axes1.set_ylabel('z (kpc)')
axes1.annotate('Plan XZ', xy=(2, 1), xytext=(-95, 92))
axes1.scatter(particles[:,0],particles[:,2],s=0.1,alpha=0.15,color="k",rasterized=True)
axes1.set_xlim(-100,100)
axes1.set_ylim(-100,100)

axes2.set_xlabel('x (kpc)')
axes2.set_ylabel('y (kpc)')
axes2.annotate('Plan XY', xy=(2, 1), xytext=(-95, 92))
axes2.scatter(particles[:,0],particles[:,1],s=0.1,alpha=0.15,color="k",rasterized=True)
axes2.set_xlim(-100,100)
axes2.set_ylim(-100,100)


fig.tight_layout()
plt.show()
print "Plot done"
plt.close()
