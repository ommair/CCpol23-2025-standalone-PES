# This is Makefile for creating executables

# This variable is for fortran compiler
FC = gfortran

# This variable will give compiler flags
FFLAGS = -O3 #-Ofast 
FCFLAGS = -c 
FLAGS = -fbounds-check

OBJS = common_arrays.o \
       vdw_ccpol23plus.o     \
       nb_induction_model.o  \
       threebody_potential.o  \
       initialization.o    \

EXECUTE = execute

$(EXECUTE): main.f90 $(OBJS)
	$(FC) $(FFLAGS) -o $(EXECUTE) main.f90 $(OBJS)
 
# These lines produces the .mod and .o files
common_arrays.o: common_arrays.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) common_arrays.f90
initialization.o: initialization.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) initialization.f90
vdw_ccpol23plus.o: vdw_ccpol23plus.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) vdw_ccpol23plus.f90 
nb_induction_model.o: nb_induction_model.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) nb_induction_model.f90
threebody_potential.o: threebody_potential.f90
	$(FC) $(FLAGS) $(FFLAGS) $(FCFLAGS) threebody_potential.f90

# This line is for clean up
clean: 
	rm -rf *.o *.mod $(OBJS) *~ $(EXECUTE) OUTPUT *.out *.xyz *.err fort.401	
