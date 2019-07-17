#!/bin/bash

source ~/.bashrc

#python3 ../testsim.py LEO0_sim_cfg.json LEO1_sim_cfg.json LEO2_sim_cfg.json LEO3_sim_cfg.json LEO4_sim_cfg.json LEO5_sim_cfg.json LEO6_sim_cfg.json sim_data.json

#python3 ../plotsim.py radar_sim_cfg.json sim_data.json


#python3 ../testodet.py LEO0_od_cfg.json LEO1_od_cfg.json LEO2_od_cfg.json LEO3_od_cfg.json LEO4_od_cfg.json LEO5_od_cfg.json LEO6_od_cfg.json sim_data.json od_output

#python3 ../plotodet.py LEO0_od_cfg.json od_output_obs_0.json od_output_estimation_0.json

#python3 ../plotodet.py LEO1_od_cfg.json od_output_obs_1.json od_output_estimation_1.json

#python3 ../plotodet.py LEO2_od_cfg.json od_output_obs_2.json od_output_estimation_2.json

#python3 ../plotodet.py LEO3_od_cfg.json od_output_obs_3.json od_output_estimation_3.json

#python3 ../plotodet.py LEO4_od_cfg.json od_output_obs_4.json od_output_estimation_4.json

#python3 ../plotodet.py LEO5_od_cfg.json od_output_obs_5.json od_output_estimation_5.json

#python3 ../plotodet.py LEO6_od_cfg.json od_output_obs_6.json od_output_estimation_6.json



#################################### ANIK SATS ####################################

#python3 ../testsim.py ANIK0_sim.json ANIK1_sim.json ANIK2_sim.json ANIK_sim_data.json

python3 ../testodet.py ANIK0_od.json ANIK1_od.json ANIK2_od.json ANIK_sim_data.json ANIK_od_output

python3 ../plotodet.py ANIK0_od.json ANIK_od_output_obs_0.json ANIK_od_output_estimation_0.json

python3 ../plotodet.py ANIK1_od.json ANIK_od_output_obs_1.json ANIK_od_output_estimation_1.json

python3 ../plotodet.py ANIK2_od.json ANIK_od_output_obs_2.json ANIK_od_output_estimation_2.json
