# OrbDetPy 

This is orbdetpy, a library of Python and Java routines for orbit determination. It is a thin Python wrapper for our estimation tools and Orekit, which are both written in Java. 

## Features

### The force model for orbit propagation currently includes:

1. EGM96 gravity field up to degree and order 360.
2. Earth solid tides due to the influence of the Sun and Moon.
3. FES 2004 ocean tide model up to degree and order 100.
4. The NRL MSISE-00 and simple exponential models for atmospheric drag.
5. Solar radiation pressure.
6. Third body perturbations from the Sun and Moon.

### The measurement model supports 
1. Range
2. RangeRate
3. Range + RangeRate
4. Azimuth + Elevation
5. RightAscension + Declination
6. PositionVelocity

Filtering is done using Orekit's Extended Kalman Filter or our custom Unscented Kalman Filter. Dynamic Model Compensation (DMC) can be used with either filter to estimate additional perturbing acclerations that result from unmodeled dynamics, maneuvers etc. **You can either use your own measurements or simulate observations using the simulateMeasurements() function**.

### The web interface currently supports

1. Select input files. Drop-down menus that allows the user to select the files in the data folder.
2. Run simulations. Runs **`testsim.py`**.
3. Add time bias, value is assigned in feature 7 or feature 10.
4. Estimate and correct observations for stellar aberration.
5. Estimate for a time bias using reference state information used for feature 5 (RA/DEC only).
6. Create a plot showing residuals between the given observations and the reference state.
7. Run orbit determination. Runs **`testodet.py`**.
8. Create plots for orbit determination analysis.
9. Run RTS UKF smoother if applicable.
10. Input desired time bias value. Input box used for feature 3.
11. Submit button to run chosen process. Clear button to reset run and allow for another run.
12. Select file to configure or view. Drop-down menu that shows files in data folder.
13. Text area that displays contents selected file and allows the user to change any of its contents (auto-save).

## Installation

OrbDetPy requires a few programs to successfully run; Python, Microsoft Visual C++ Build Tools, and JDK 8+.

1. Install [Python 3.6+](https://www.python.org/downloads/). To see if you have already installed Python 3.6+ or to see if it installed correctly, go to command prompt and enter `python`. Your current working version of python should be displayed. If any errors occurred, download/re-install python from the link provided and make sure to select ADD TO PATH during the installation procedure. 
2. Install `pyjnius`. This requires Microsoft Visual C++ Build Tool, Java Development Kit 8+, and `Cython` to install. An in-depth installation guide for `pyjnius` can be found [here](https://pyjnius.readthedocs.io/en/stable/installation.html).
3. Now pip install the following packages 
```bash
pip install numpy scipy pyjnius matplotlib requests astropy plotly dash
```
4. Finally, clone the orbdetpy repo to your local computer. For instruction on how to clone repositories, follow this [guide](https://help.github.com/en/articles/cloning-a-repository).
5. For those still having issues, ensure you are following the installation procedure in order and exactly as specified. Also ensure installs are all 32-bit or all 64-bit.

The orbdetpy/lib folder contains JAR files for the following libraries, which are imported by orbdetpy automatically.

1. [Google GSON](https://github.com/google/gson)
2. [Hipparchus](https://hipparchus.org)
3. [Orekit](https://www.orekit.org)

The astrodynamics data in orbdetpy/data must be updated periodically using examples/update_data.py.

## Command Line Examples

The following example programs can be found in the 'examples' folder. These examples use the Python wrapper interface but calling the underlying Java implementation directly is straightforward.

1. **`testsim.py`** : Demonstrates the measurement simulator. Note that
   maneuvers can be incorporated into the force model during simulation.
2. **`plotsim.py`** : Plots the results of simulations created using `testsim.py`.
3. **`testodet.py`**: Demonstrates orbit determination in orbdetpy.
4. **``plotodet.py``** : Plots the results of fitting orbits using `testodet.py`.
5. **`run_tests.py`** : Run all the use cases under examples/data. Simulated measurements, orbit fits, differences between simulated truth versus estimates, and 3-sigma of estimated covariances will be written to output/ sub-folders.

orbdetpy uses JSON files to store settings, measurements and estimation results. The files in examples/data show how to configure measurement simulation and orbit determination using radar or telescope data. The file docs/file_formats.md documents the structure of the JSON files. The following are some typical use cases. It is assumed that the current working directory is examples/data.

1. Simulate state vectors and radar measurements
```bash
   python ../testsim.py radar_sim_cfg.json sim_data.json
```
This will run the simulation configured in radar_sim_cfg.json and write simulated output to sim_data.json.

2. Plot simulation results:
```bash
   python ../plotsim.py radar_sim_cfg.json sim_data.json
```
This will plot the simulated data generated in (1).

3. Run OD on simulated radar data:
```bash
   python ../testodet.py radar_od_cfg.json sim_data.json od_output.json
```
This will run OD on the simulated radar data generated in (1) using the OD configuration in radar_od_cfg.json and write OD output to od_output.json.

4. Plot OD results:
```bash
   python ../plotodet.py radar_od_cfg.json sim_data.json od_output.json
```
   This will plot the OD results from (3).

**To learn more about the configure file format, go to `orbdetpy/docs/filefmt.rst`**

## Web Interface Examples

An alternative to the above example is to use the web interface. Executing ``runWebInterface.py`` will create a locally hosted site. The web address will display in the terminal. Open a browser of your choice and navigate to the displayed address. Follow the instructions in the About section, or follow the instructions below for a simple use case.

1. Simulate state vectors and radar measurements
    1. Select the "Select Simulation Configure File" dropdown: `radar_sim_cfg.json`
    2. Select the "Run Simulation" checkbox:
    3. Press the "Submit" button:
    4. Once complete, press the "Clear" button: This will reset the program and allow for another process to execute.
   
This will run the simulation configured in `radar_sim_cfg.json` and write simulated output to `radar_obs_data.json`.
   
2. Run OD on simulated radar data

    1. Select the "Select OD Configure File" and "Select Observation File" dropdown: radar_od_cfg.json, radar_obs_data.json.
    2. Select the "Run orbit determination" checkbox:
    3. Press the "Submit" button:
    4. Once complete, press the "Clear" button: This will reset the program and allow for another process to execute.

This will run OD on the simulated radar data generated in (1) using the OD configuration in `radar_od_cfg.json` and write OD output to `radar_od_output.json`.

3. Plot results ( default save directory `orbdetpy/examples/Plots` ).

    1. Select the "Select OD Configure File", "Select Observation File", and "Select OD Output File" dropdown: `radar_od_cfg.json`, `radar_obs_data.json`, `radar_od_out.json`.
    2. Select the "Create plots for orbit determination analysis" checkbox.
   3. Press the "Submit" button.
   4. Once complete, press the "Clear" button: This will reset the program and allow for another process to execute.

This will plot OD results from (2).

**To learn more about the configure file format, go to `orbdetpy/docs/filefmt.rst`**

## Future Work

The following tasks are under consideration. Community contributions are
always welcome.

1. A batch least squares implementation.
2. Parametric analysis i.e. the ability to pass-through certain
   measurement types.
3. Time bias estimation during filter run.
4. Dynamic outlier filtering

## Bug Reports

Comments, criticisms and bug reports are very welcome and may be sent to the package maintainer by email or the project's website.

Shiva Iyer <shiva.iyer AT utexas DOT edu>
Additional Contact:
Vivek Desai <vivek.desai AT utexas DOT edu>