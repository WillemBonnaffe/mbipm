# mbipm

## Get started

First make sure that you have installed Python and Anaconda for managing virtual environments.

To install the environment required for running the model run the following commands:

```
conda env create --name mbipm --file=environments.yml
conda activate mbipm
```

If this does not work you may have to create an empty environment and install packages manually:

```
condra create env --name mbipm python=3.10
conda activate mbipm
conda install pytorch
conda install numpy
conda install matplotlib
conda install seaborn
conda install 
```
