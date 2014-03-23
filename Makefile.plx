###############################################################################
# Rules
#
all: pasta_axi.exe
default: pasta_axi.exe

include Dependencies.inc
# Is Dependencies.inc available ?
Dependencies.inc:
	@echo "##################################################################"
	@echo BEFORE COMPILING, YOU SHOULD HAVE AN APPROPRIATE FILE
	@echo Dependencies.inc AVALAIBLE.
	@echo "##################################################################"
	@exit 1

release:
	@echo making release
	@rm -f *.mod *.o
	@cd .. && tar -czf ./backups/pasta_axi_$(NOW).tar.gz pasta_axi/

pasta_axi.exe: main.o $(OBJ) 
	$(F90LNK) $(F90FLAGS) main.o $(OBJ) $(LIBS) -o $(FOLDER)$(NAME)

clean:
	@echo cleaning
	@rm -f *.o *.mod

###############################################################################
# Define variables

NOW := $(shell date +"%c" | tr ' :' '__')

# Compilers
FC = mpif90
CC = mpicc

F90CMP = $(FC) -c
F90LNK = $(FC)
F90OPT = -O3 -fopenmp
F90DEB = #-g -Wall #-fcheck=all
F90FLAGS = $(F90OPT) $(F90DEB)

CCMP   = $(CC) -c
CLNK   = $(CC)
COPT   = -O3
CDEB   =
CFLAGS = $(COPT) $(CDEB)

# Include paths
F90INCDIR = -I$(OPENMPI_HOME)/include -I$(MUMPS_HOME)/include -I$(MKL_INC)
CINCDIR   = -I libLOCA

# Libraries
LIBARPACK = -L/cineca/prod/libraries/arpack/96/gnu--4.5.2/lib -larpack
LIBLAPACK = -L$(LAPACK_LIB) -llapack
LIBBLAS   = -L$(BLAS_LIB) -lblas
LIBSCOTCH = -L$(SCOTCH_HOME)/lib -lptscotch -lptesmumps -lptscotcherr
LIBMKL    = -L${MKL_LIB} -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64
LIBMUMPS  = -L$(MUMPS_LIB) -ldmumps -lzmumps -lpord -lmumps_common -lpord
LIBMPI    = -L$(OPENMPI_HOME)/lib -lmpi
LIBLOCA   = -L./libLOCA/ -lloca
LIBS      = $(LIBMUMPS) $(LIBMPI) $(LIBLOCA) $(LIBARPACK) $(LIBSCOTCH) $(LIBMKL) $(LIBLAPACK) $(LIBBLAS) $(LIBMUMPS) $(LIBSCOTCH)

FOLDER = $(HOME)/software/bin/
NAME   = pasta_axi


###############################################################################
# Rules

.SUFFIXES:
.SUFFIXES: .c .f90 .o

.f90.o:
	$(F90CMP) $(F90FLAGS) $(F90INCDIR) $<

.c.o:
	$(CCMP) $(CFLAGS) $(CINCDIR) $<
