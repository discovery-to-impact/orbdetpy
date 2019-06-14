#!/bin/bash

source ~/.bashrc

#python3 ../testsim.py radec_sim_cfg.json sim_data.json

#python3 ../plotsim.py radar_sim_cfg.json sim_data.json

python3 ../testodet.py radec_od_cfg.json sim_data.json od_output.json

python3 ../plotodet.py radec_od_cfg.json sim_data.json od_output.json

