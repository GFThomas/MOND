 **** EXTRACT AND ANALYSE RAMSES/POR DATA ****

Download instructions:
git clone https://github.com/GFThomas/MOND.git
Please then cut the extract_por folder and paste it into another directory of your choice. You may then delete the MOND folder if you wish.

Files inside TAR archive (algorithm consists of the following 17 files):

README
Makefile
rnd_array90.f90
xybinsuite.f
deviates.f90
rqsort.f90
densitycenter90.f90
commons4ramses.f90
format4ramses.f90
radial_profile.f90
partial_com.f90
xyimage.f90
vec3routines.f
xbin05.f
nextline.f
analyse_ramses_data.f90
fmtramses.par

Installation:
make xpordata

Running:
./xpordata.x

The fmtramses.par needs to be in the same folder from which the user calls the above command unless specified as an argument:
./xpordata.x <yourparamsfile>

The results are also written in this folder.

Basic options in the parameter file which may be frequently changed (line number of parameter file given at start):
1  Relative path to RAMSES simulation output folder (use '.' if same folder)
2  Output No. (negative for all up to there, e.g. '-45' means all from 1 to 45). Using e.g. 45 will extract just snapshot 45, or however many snapshots there are (if this is <45).
3  Number of CPU threads used by RAMSES during the simulation (leave at <= 0 for automatic detection)
...(Other options described in parameters file)
8  Scale position unit (PoR output is multiplied by this factor)
9  Scale velocity unit (as above)
10 Scale mass unit (as above)
11 Scale acceleration (if exists) and birth time (if exists) units.
12 Boolean flag for existence of accelerations block (it exists in the 2015 POR
   version but not in standard RAMSES).
13 Optional rotation of the ASCII output (use 0. 0. 0. for no rotation)
...
15 first number: Select mode: 0 - only convert to ASCII
                1 - only analyse (assumes ASCII exists)
		>1- do both
	 	<0- do neither (e.g. just create gnuplot image data from existing ASCII)
If only extracting data to ASCII, options beyond this point are unused and should be left as they are.
15 second number: Volume dimension used for binning (e.g. surface or spatial density, 2 = cylindrical binning, 3 = spherical binning)

The remaining options are for radial binning (no warranty!), image (gnuplot format),
and for a special kind of dark matter simulations (non-MONDian RAMSES) with three halos
used in Oehm, Thies & Kroupa (2017, DOI: 10.1093/mnras/stw3381).
For general purpose use they may be ignored safely.
The only important thing might be to type the box length into line 20 (first entry),
then data will be relative to box centre for binning only (not ASCII output above).
Lines below this and lines not described here should be left unaltered.
Line ~190 in format4ramses.f90 is the really important line for normal options
(write ASCII files without acceleration and birth time). It is marked with a similar comment.
