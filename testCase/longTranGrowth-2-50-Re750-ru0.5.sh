#!/bin/bash
#
# simple script to automate the computation of the transient energy growth
# for various taus. This script has to be used because reading taus from
# 'taus.in' doesn't work yet (there's a bug somewhere).
#
# author: jacopo.canton@mail.polimi.it

mkdir tranGrowthOut/ru0.5/Re750

Time="2"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="5"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="8"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="10"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="20"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="30"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="40"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

# sleep 30

Time="50"
sed -i "/tau/c\\$Time     # tau (end time)" program_data.in
mpirun -np 4 pasta_axi > log`echo "($Time)" | bc -l`

mv log* tranGrowthOut/ru0.5/Re750/
mv noh* tranGrowthOut/ru0.5/Re750/
mv tranGrowthOut/*.vtk tranGrowthOut/ru0.5/Re750/
mv tranGrowthOut/*.dat tranGrowthOut/ru0.5/Re750/
