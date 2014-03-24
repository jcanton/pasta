############################################################
#
# Compiler and linker
#

CMP = gfortran -c
LNK = gfortran
OPT = -O2


############################################################
#
# Objects: list of all objects *.o

OBJS =

############################################################
# Executable generation
#

coaxialTubes:     coaxialTubes.o  $(OBJS)
	$(LNK) $(OPT)  coaxialTubes.o  $(OBJS) \
                 -o coaxialTubesBuildMesh

############################################################
# Objects generation

coaxialTubes.o:   coaxialTubes.f90  $(OBJS)
	$(CMP) $(OPT) coaxialTubes.f90

############################################################
# Cleaning command to remove all objects *.o, files *.mod and
# the executables *.exe

clean: 
	@echo cleaning
	@rm -f  *.o  *.mod coaxialTubesBuildMesh
  
# The command "make clean" deletes all the indicated files

############################################################
