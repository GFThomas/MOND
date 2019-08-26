F90=gfortran
F77=gfortran
CC=gcc

# standard fast flags
#FFLAGS= -O2
#FFLAGS= -fbounds-check
FFLAGS= -fbounds-check -O2
#FFLAGS= -g -pedantic -Wall -fbounds-check
#FFLAGS= -pedantic -Wall -fbounds-check

DDFLAGS= -fno-range-check -ffree-form

XPORDATA=rnd_array90.o xybinsuite.o deviates.o rqsort.o densitycenter90.o commons4ramses.o format4ramses.o radial_profile.o partial_com.o xyimage.o vec3routines.o xbin05.o nextline.o analyse_ramses_data.o

.SUFFIXES: .f90 .F90 .f .F .f77 .F77
.f90.o:
	$(F90) $(FFLAGS) -c $<
.F90.o:
	$(F90) $(FFLAGS) $(CFLAGS) -c $<
.f.o:
	$(F77) $(FFLAGS) -c $<
.F.o:
	$(F77) $(FFLAGS) -c $<
.f77.o:
	$(F77) $(FFLAGS) -c $<
.F77.o:
	$(F77) $(FFLAGS) -c $<

xpordata: $(XPORDATA)
	$(F90) $(FFLAGS) -o xpordata.x $(XPORDATA)

clean:
	\rm *.o *.mod

purge:
	\rm *.o *.mod *.x

rmtilde:
	\rm *~

ppurge:
	\rm *.o *.mod *.x \rm *~
