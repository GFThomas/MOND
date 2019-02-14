# Pympor


## About

Python utils for read POR outputs

Fev. 2015 by Guillaume THOMAS <guillaume.thomas@astro.unistra.fr>
Copyright (c) 2015 Observatoire de Strasbourg. All rights reserved.




### Packages needed

* `python`
* `gfortran` or other fortran compiler that works with f2py
* `f2py` <https://sysbio.ioc.ee/projects/f2py2e/>

**Python packages**
* `numpy` <http://www.numpy.org/>
* `mathplotlib`<http://matplotlib.org/>




### Installation :
	
	>make
	>echo "export PYTHONPATH=\$PYTHONPATH:$(pwd)" >> ~/.bashrc




### Examples :

This package contain 2 examples, with the outputs comming from one of the simulation from Thomas et al. 2018a

	> cd example


** The first example is to plot the particle from the 2nd output

	>  python plot_particles.py
		2
		
		
** The second is to plot the level of refinement of the AMR grid around the globular cluster on the plan XZ (can take a long time)

	> python plot_amr.py
		2
		1
		1
		7.0
		-30
		30

The choice of 7.0 along the y-axis is for speed reason, if you want to see the full AMR around globular cluster, you should try with 7.9





### Programs in that package

`Pympor.py` Program for read POR outputs.

	* chargement_part : Read the particles of the POR or Ramses outputs
	                    # INPUT : n° of the POR or Ramses output
	                    # OUTPUT : position(x,y,z in kpc), velocity(vx,vy,vz in km/s), mass(M☉) and n° of each particles and the time of the output

	* chargement_grav : Read the grav file of the POR or Ramses outputs
	                    # INPUT : n° of the POR or Ramses output and 1 if MOND and 0 if Newton
	                    # OUTPUT : position(x,y,z in kpc), size (kpc), potential(m²/s²), acceleration(ax,ay,az in m/s²), AMR level of each grid, newtonian potential(m²/s²), newtonian acceleration(ax,ay,az in m/s²)

	* filters : Select particles abrove and below a mass m_threshold
	                    # INPUT : particles, m_threshold
	                    # OUTPUT : particles abrove and particles below m_threshold


 

`data/dumpdmparts.f90` fortran program to read the part files


`data/gravtoascii.f90` fortran program to read the grav files
If you have some problem with your RAM memory you can decrease the xxxxxxx value in line 62 :
  !=============================================================
  character(LEN=12),intent(in)::name_file
  integer,intent(in):: log_mond
  integer::id_grav_out
  real(KIND=8),dimension(1:xxxxxxx,1:9),intent(out)::grav_out
  !=============================================================
