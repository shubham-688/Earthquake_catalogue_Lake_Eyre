#
#	Make file for program shake_na.f
#

# location of fortran compiler

#F77 = /opt/local/SUNWspro/bin/f77

# compiler options

FFLAGS = -C

#---- start user specific make --------------------------------------------
# user supplied subroutines called by user interface routines

SDR = shaklib
SHK = $(SDR)/libsun.o $(SDR)/libtau.o \
      $(SDR)/libshaknd.o $(SDR)/libshaki.o $(SDR)/ellip.o

#------------------------------------------------------------------------

all: na_subs shake_na 


shake_na:
	gfortran $(FFLAGS) -o ../bin/shake_na shake_na.f \
	lib/na.o $(SHK)
	cp shake_na.f Backups 

na_subs:
	gfortran $(FFLAGS) -c na.f
	mv na.o lib

shake_msf:
	gfortran $(FFLAGS) -o ../bin/shake_msf shake_msf.f \
	$(SHK)

#------------------------------------------------------------------------
