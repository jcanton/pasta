#!/bin/bash
#
# make utility for PaStA
# jcanton@mech.kth.se

###############################################################################
# define local settings
#
srcDir="src/"
makefile="Makefile.erebos"
#
###############################################################################
# set parameters and make
#
export CASEDIR="`pwd`"
export PASTANAME="pasta"
export PASTAFLAGS="" #-DASCIIEIGENVECTOR -DASCIIRESTART -DSAVEEIGENVECTOR

cd ${srcDir}
make -f ${makefile}
