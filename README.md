# MBIPM

This repository contains the code necessary to build, simulate, and fit multi-species bioenergetic integral projection models.

## Get started

First make sure that you have installed Python and Anaconda for managing virtual environments.

To install the environment required for running the model run the following commands:

```
conda env create --name mbipm --file=environment.yml
conda activate mbipm
```

If this does not work you may have to create an empty environment and install packages manually:

```
conda create -n mbipm python=3.10 numpy matplotlib seaborn jupyter
conda activate mbipm
pip3 install torch torchvision
```
Now you can open jupyter lab by running in the terminal the following command:

```
jupyter lab
```

You are ready to explore the notebooks.


## Notebooks

The repository contains multiple notebooks that can be run sequentially.

`A1-format-time-series.r` is an R script that formats the full time series of animal counts in Yellowstone.

`C1-setup-calibration.ipynb` is a jupyter notebook that prepares the Yellowstone ecosystem model (YEM) for fitting to the time series.  

`C2-calibrate.ipynb` fits the YEM to the time series using differential-evolution (DE) optimisation.

`C3-check-calibration.ipynb` checks the quality of the fit of the model to the time series.

`D1-simulate.ipynb` takes the calibrated model and performs long term simulations.

`D2-simulate-extirpation.ipynb` takes the calibrated model and performs long term simulation following predator extirpation and reintroduction.


## Supporting files

The notebooks call various files which are detailed below.

`data/master-v2025-05-15-17-06.xlsx` contains all raw time series of elk, bison, cougar, and wolves.

`À1-outputs/A1-tab-all-time-series.csv` is the formatted time series that can be used to calibrate the model.

`f_YEM_<model-version-name>.py` defines all functions necessary to initiate and simulate the YEM.

`f_YEM_wrapper_<model-version-name>.py` defines a wrapper to optimise the paramters of the YEM following DE optimisation. 

`f_demc_sparse_noise.py` contains functions that are necessary to perform DE optimisation.

`ìnitial-conditions-<model-version-name>.yaml` contains the initial state of the variables in the YEM before calibration.

`ìnitial-conditions-apriori-<model-version-name>.yaml` contains the initial state of the variables in the YEM after calibration.

`parameters-aposteriori-<model-version-name>.yaml` contains the values of all parameters in the YEM after calibration.

`parameters-apriori-<model-version-name>.yaml` contains the values of all parameters in the YEM before calibration.

`settings-apriori-<model-version-name>.yaml` contains the settings (e.g. number of bins) of the YEM both for usable before and after calibration.

