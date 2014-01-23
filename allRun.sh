#!/bin/bash
#
# simple script to automate the computation of the transient energy growth
# for various taus. This script has to be used because reading taus from
# 'taus.in' doesn't work yet (there's a bug somewhere).
#
# author: jacopo.canton@mail.polimi.it

t0="6"
tf="26"
dt="2"

Time=$t0

sed -i "/tau/c\\$Time     # tau (end time)" program_data.in

while [ $Time -lt $tf ]; do

   Time=`expr $Time + $dt`

   mpirun -np 4 continuation_axi_par > log`echo "$Time-$dt"` && sed -i "/tau/c\\$Time # tau (end time)" program_data.in

done
