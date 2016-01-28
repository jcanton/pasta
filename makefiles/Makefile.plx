###############################################################################
module purge
module load profile/advanced
module load autoload mumps
module load arpack/96--gnu--4.5.2
module load lapack/3.3.1--gnu--4.5.2
module load blas/2007--gnu--4.5.2
module load scalapack/2.0.2--openmpi--1.6.3--gnu--4.5.2
###############################################################################
# Rules
#
all: pasta_loca.exe
default: pasta_loca.exe

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
	@cd .. && tar -czf ./backups/pasta_loca_$(NOW).tar.gz pasta_loca/

pasta_loca.exe: main.o $(OBJ) 
	$(F90LNK) $(F90FLAGS) main.o $(OBJ) $(LIBS) -o $(FOLDER)$(NAME)

clean:
	@echo cleaning
	@rm -f *.o *.mod

###############################################################################
# Define variables

NOW := $(shell date +"%c" | tr ' :' '__')

# Compilers
FC = mpif90 -cpp

F90CMP = $(FC) -c
F90LNK = $(FC)
F90OPT = -O3 -fopenmp
F90DEB = -DASCIIEIGENVECTOR -DDEBUG=0 -DMUMPSDEBUG=0 #-g -Wall #-fcheck=all -DASCIIRESTART
F90FLAGS = $(F90OPT) $(F90DEB)

# Include paths
F90INCDIR = -I$(OPENMPI_HOME)/include -I$(MUMPS_HOME)/include -I$(MKL_INC)

# Libraries
LIBARPACK = -L$(ARPACK_LIB) -larpack
LIBLAPACK = -L$(LAPACK_LIB) -llapack
LIBSCALAPACK = -L$(SCALAPACK_LIB) -lscalapack
LIBBLAS   = -L$(BLAS_LIB) -lblas
LIBSCOTCH = -L$(SCOTCH_HOME)/lib -lptscotch -lptesmumps -lptscotcherr
LIBMKL       = -L$(MKL_LIB) -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm
LIBMUMPS  = -L$(MUMPS_LIB) -ldmumps -lzmumps -lpord -lmumps_common -lpord
LIBMPI    = -L$(OPENMPI_HOME)/lib -lmpi
LIBS      = $(LIBMUMPS) $(LIBMPI) $(LIBARPACK) $(LIBSCOTCH) $(LIBSCALAPACK) $(LIBMKL) $(LIBLAPACK) $(LIBBLAS) $(LIBMUMPS) $(LIBSCOTCH)

FOLDER = ./
NAME   = pasta_loca

###############################################################################
# Rules

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(F90CMP) $(F90FLAGS) $(F90INCDIR) $<
