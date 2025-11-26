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
conda create -n mbipm python=3.10 numpy matplotlib seaborn jupyter
conda activate mbipm
conda install pytorch cpuonly -c pytorch 
```
Now you can open jupyter lab by running in the terminal the following command:

```
jupyter lab
```

You are ready to explore the notebooks.

