=================================================
orbdetpy - Orbit determination in Python and Java
=================================================

Introduction
------------

This is orbdetpy, a library of Python and Java routines for orbit
determination. It is a thin Python wrapper for our estimation tools
and Orekit, which are both written in Java. 

Features
--------

The force model for orbit propagation currently includes:

1) EGM96 gravity field up to degree and order 360.
2) Earth solid tides due to the influence of the Sun and Moon.
3) FES 2004 ocean tide model up to degree and order 100.
4) The NRL MSISE-00 and simple exponential models for atmospheric drag.
5) Solar radiation pressure.
6) Third body perturbations from the Sun and Moon.

The measurement model supports range, range-rate, angles, and inertial
Cartesian coordinates. Filtering is done using Orekit's Extended Kalman
Filter or our custom Unscented Kalman Filter. Dynamic Model Compensation
(DMC) can be used with either filter to estimate additional perturbing
acclerations that result from unmodeled dynamics, maneuvers etc.

You can either use your own measurements or simulate observations using
the simulateMeasurements() function.

The web interface for process execution currently includes:

About – provides brief information on the interface
1) Select input files. Multi option dropdown that allows the user to select the files in the data folder.
2) Data segmentation. Range slider that allows the user to select the time range of data to run. Displays the start and end time below the slider. (in-progress)
3) Run simulation. Runs testsim.py
4) Parse data. Parses currently provided L3 ADS data
5) Add time bias, value is assigned in feature 7 or feature 12
6) Estimate and correct observations for stellar aberration
7) Estimate for a time bias using reference state information used for feature 5 (only works for RA/DEC right now)
8) Create a plot showing residuals between the given observations and the reference state
9) Run orbit determination. Runs testodet.py
10) Run RTS UKF smoother (must run with feature 11 for now)
11) Create plots for orbit determination analysis
12) Input desired time bias value. Input box used for feature 5
13) Submit button to run chosen process. Clear button to reset run and allow for another run.
14) Select file to configure or view. Dropdown that shows files in data folder
15) Text area that displays contents selected file and allows the user to change any of its contents

Installation
------------

1) Python 3.6+ must be installed with the packages numpy, scipy, pyjnius,
   matplotlib, requests, plotly, and dash.
2) Install the Java Development Kit 8+ (1.8+) from `here
   <http://openjdk.java.net>`_. Set the JAVA_HOME environment variable
   to point to your JDK installation.

The lib/ folder contains JAR files for the following libraries, which are
imported by orbdetpy automatically.

1) `Google gson <https://github.com/google/gson>`_
2) `Hipparchus 1.4+ <https://hipparchus.org>`_ 
3) `Orekit 9.3+ <https://www.orekit.org>`_

Space weather data in data/ is obtained from `CelesTrak <http://www.celestrak.com/SpaceData/>`_.
The astrodynamics data in data/ must be periodically updated using the
script examples/update_data.py.

Command Line Examples
---------------------

The following example programs can be found in the 'examples' folder.
These examples use the Python wrapper interface but calling the
underlying Java implementation directly is quite straightforward.

1) testsim.py : Demonstrates the measurement simulator. Note that
   maneuvers can be incorporated into the force model during simulation.

2) plotsim.py : Plots the results of simulations created using testsim.py.

3) testodet.py : Demonstrates orbit determination in orbdetpy.

4) plotodet.py : Plots the results of fitting orbits using testodet.py.

5) runWebInterface : Executes locally hosted web interface. 
orbdetpy uses JSON files to store settings, measurements and estimation
results. The files in examples/data show how to configure measurement
simulation and orbit determination using radar or telescope data. The
file docs/filefmt.rst documents the structure of the JSON files.

The following are some typical use cases. It is assumed that the current
working directory is examples/data.

1) Simulate state vectors and radar measurements:

   python ../testsim.py radar_sim_cfg.json sim_data.json

   This will run the simulation configured in radar_sim_cfg.json and
   write simulated output to sim_data.json.

2) Plot simulation results:

   python ../plotsim.py radar_sim_cfg.json sim_data.json

   This will plot the simulated data generated in (1).

3) Run OD on simulated radar data:

   python ../testodet.py radar_od_cfg.json sim_data.json od_output.json

   This will run OD on the simulated radar data generated in (1)
   using the OD configuration in radar_od_cfg.json and write OD
   output to od_output.json.

4) Plot OD results:

   python ../plotodet.py radar_od_cfg.json sim_data.json od_output.json

   This will plot the OD results from (3).

Web Interface Examples
----------------------

An alternative to the above example is to use the web interface. Executing 
runWebInterface.py will create a locally hosted site. The web address will 
display in the terminal. Open a browser of your choice and navigate to the
displayed address. Follow the instructions in the About section, or 
follow the instructions below for a simple use case.

1) Simulate state vectors and radar measurements:
   
   Select the "Select input files" dropdown: radar_sim_cfg.json
   Select the "Run Simulation" checkbox:
   Press the "Submit" button:
   This will run the simulation configured in radar_sim_cfg.json and
   write simulated output to radar_obs_data.json.
   Once complete, press the "Clear" button:
   This will reset the program and allow for another process to execute.
   
2) Run OD on simulated radar data:

   Select the "Select input files" dropdown: radar_od_cfg.json, radar_obs_data.json
   Select the "Run orbit determination" checkbox:
   Press the "Submit" button:
   This will run OD on the simulated radar data generated in (1)
   using the OD configuration in radar_od_cfg.json and write OD
   output to od_output.json.
   Once complete, press the "Clear" button:
   This will reset the program and allow for another process to execute.

3) Plot results:

   Select the "Select input files" dropdown: radar_od_cfg.json, radar_obs_data.json, radar_od_out.json
   Select the "Create plots for orbit determination analysis" checkbox:
   Press the "Submit" button:
   This will plot OD results from (2).
   Once complete, press the "Clear" button:
   This will reset the program and allow for another process to execute.

Future Work
-----------

The following tasks are under consideration. Community contributions are
always welcome.

1) A batch least squares implementation.
2) Parametric analysis i.e. the ability to pass-through certain
   measurement types.
3) Time bias estimation during filter run.
4) Dynamic outlier filtering

Bug Reports
-----------

Comments, criticisms and bug reports are very welcome and may be sent to
the package maintainer by email or the project's website.

Shiva Iyer <shiva.iyer AT utexas DOT edu>
