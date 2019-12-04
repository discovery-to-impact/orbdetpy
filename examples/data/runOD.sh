#!/bin/bash

source ~/.bashrc



#python3 ../testsim.py ANIK0_sim.json ANIK1_sim.json ANIK2_sim.json LEO0_sim.json LEO1_sim.json sim_data.json

#python3 ../testodet.py ANIK0_od.json ANIK1_od.json ANIK2_od.json LEO0_od.json LEO1_od.json sim_data.json od_output

#python3 ../testsim.py MEO0_sim.json MEO1_sim.json MEO2_sim.json sim_data.json

#python3 ../testodet.py MEO0_od.json sim_data.json od_output




python3 ../plotodet.py MEO0_od.json od_output_obs_0.json od_output_estimation_0.json

python3 ../plotodet.py MEO0_od.json od_output_obs_1.json od_output_estimation_1.json

python3 ../plotodet.py MEO0_od.json od_output_obs_2.json od_output_estimation_2.json
