my_objs     = externals_all.o get_cc_info.o

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
