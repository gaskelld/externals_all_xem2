my_objs     = externals_all.o readin_ext.o radiators.o continuum.o atailfl.o qetail.o qepeak.o epeak.o func.o cross.o elastic.o qelastic.o quasiy8.o fyxsec8.o secnuclw.o inelast.o ineft.o f2glob.o f2allm.o f1f207.o f1f2qe07.o f1f209.o f1f2qe09.o gsmearing.o f1f2qe16.o inelastic_models.o get_cc_info.o

#CERNLIBS = /apps/cernlib/x86_64_rhel7/2005/lib -lmathlib 
#CERNLIBS = /u/site/cernlib/x86_64_rhel7/2005/lib -lmathlib
#CERNLIBS = /opt/cernlib/2005/lib -lmathlib
CERNLIBS = $(CERN_ROOT)/lib -lmathlib

FFLAGSA = -C -g -w -fno-automatic -fbounds-check -ffixed-line-length-132 -finit-local-zero 
INCLUDES = -I.
FFLAGS = $(INCLUDES) $(FFLAGSA)
FFLAG1 = $(FFLAGS) -c
OTHERLIBS = \
-L$(CERNLIBS) 
FC  := gfortran
CPP := g++

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

all: externals_all

externals_all :  Makefile $(my_objs)
	$(FC) -o $@ $(FFLAGS) $(my_objs)  $(OTHERLIBS)
	cp externals_all ../.
clean:
	rm -f *.o *.so externals_all ../externals_all
