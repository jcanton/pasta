#!/bin/bash
#
# simple script to automate the computation of the transient energy growth
# for various taus. This script has to be used because reading taus from
# 'taus.in' doesn't work yet (there's a bug somewhere).
#
# author: jacopo.canton@mail.polimi.it


# Time="2"
# mkdir tranGrowthOut/Re750/evolveTau2
# sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
# sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau2\/\' # dns_output_directory" program_data.in
# sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_002000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
# mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &
# 
# sleep 30
# 
# Time="5"
# mkdir tranGrowthOut/Re750/evolveTau5
# sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
# sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau5\/\' # dns_output_directory" program_data.in
# sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_005000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
# mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &
# 
# sleep 30
# 
# Time="8"
# mkdir tranGrowthOut/Re750/evolveTau8
# sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
# sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau8\/\' # dns_output_directory" program_data.in
# sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_008000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
# mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &
# 
# sleep 30
# 
# Time="10"
# mkdir tranGrowthOut/Re750/evolveTau10
# sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
# sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau10\/\' # dns_output_directory" program_data.in
# sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_010000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
# mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &
# 
# sleep 30

Time="20"
mkdir tranGrowthOut/Re750/evolveTau20
sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau20\/\' # dns_output_directory" program_data.in
sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_020000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &

sleep 30

Time="30"
mkdir tranGrowthOut/Re750/evolveTau30
sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau30\/\' # dns_output_directory" program_data.in
sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_030000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &

sleep 30

Time="40"
mkdir tranGrowthOut/Re750/evolveTau40
sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau40\/\' # dns_output_directory" program_data.in
sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_040000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l` &

sleep 30

Time="50"
mkdir tranGrowthOut/Re750/evolveTau50
sed -i "/dns_end_time/c\\$Time # dns_end_time" program_data.in
sed -i "/dns_output_directory/c\\'tranGrowthOut\/Re750\/evolveTau50\/\' # dns_output_directory" program_data.in
sed -i "/optimal/c\\'tranGrowthOut\/Re750\/tranGrowthShape0-tau_050000SteadyState\.vtk\' # optimal initial perturbation" program_data.in
mpirun -np 4 pasta_axi > evolveLog`echo "($Time)" | bc -l`
