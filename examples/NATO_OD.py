# run_tests.py - Program to run simulation and OD tests.
# Copyright (C) 2019 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import json
import math
import time
from threading import Semaphore
from orbdetpy.estimation import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot

def callback(data, extra):
    if (isinstance(data, Exception)):
        print(data)
    sem.release()

print("run_tests start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
odpdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
datdir = os.path.join(odpdir, "examples", "data")

lock_count = 0
sem = Semaphore(0)
for root, dirs, files in os.walk(datdir):
    outdir = os.path.join(root, "output")
    if ("output" in dirs):
        dirs.remove("output")

    for fname in files:
        idx = fname.find("od_cfg.json")
        if (idx == -1):
            continue

        print("Fitting {}".format(fname))
        with open(os.path.join(outdir, fname[:idx] + "obs.json"), "r") as fp:
            obs = json.load(fp)

        for algo in ["EKF","UKF"]:
            with open(os.path.join(root, fname), "r") as fp:
                config = json.load(fp)
                config["Estimation"]["Filter"] = algo
            lock_count += 1
            determine_orbit(config, obs, async_callback = callback,
                            output_file = os.path.join(outdir, fname[:idx] + algo + "_fit.json"))

for i in range(lock_count):
    sem.acquire()

plot("/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/1999-055A_od_cfg.json",
        "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/1999-055A_obs.json",
        "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/1999-055A_EKF_fit.json"
        ,interactive = True)

plot("/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/1999-055A_od_cfg.json",
        "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/1999-055A_obs.json",
        "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/1999-055A_UKF_fit.json"
        ,interactive = True)

print("run_tests end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
