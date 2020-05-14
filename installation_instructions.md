# Installation instructions

The Santa Cruz County COVID-19 model is written in [Stan](https://mc-stan.org/). It comes with a [Jupyter notebook](https://jupyter.org/) written in Python that performs some data processing and visualization. The instructions below provide basic guidance for installing Python, Jupyter notebooks, and [pystan](https://mc-stan.org/users/interfaces/pystan) (the Python interface to Stan). 

*NOTE:* The Stan code for the model can also be run in R (using [RStan](https://mc-stan.org/users/interfaces/rstan)) but we do not provide a Jupyter notebook written in R.

## Installation instructions for macOS

We recommend to install Anaconda via the instructions [here](https://docs.anaconda.com/anaconda/install/). The Anaconda installation will provide Jupyter notebooks along with Python, including the `numpy` and `scipy` packages.

### pystan 

Once Anaconda is installed, run the follwing shell command to install pystan (see [here](https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/) for more information):
```
conda install pystan 
```


## Installation instructions for Linux

If it is not already installed, install Python 3 and pip. In Ubuntu, use the shell command:
```
sudo apt install python3 python3-pip
```
[This article](https://www.tecmint.com/install-pip-in-linux/) describes the installation of pip for a variety of Linux distributions. If a Python 3 installation from source is required (typically it is not), instructions can be found [here](https://solarianprogrammer.com/2017/06/30/building-python-ubuntu-wsl-debian/)

Now, the required python packages can be installed using pip, the Python package installer. In Ubuntu, `pip3` is the Python 3 version of pip. 

### Jupyter notebooks

After installing, run the following shell command to install jupyter notebooks (you may need to use `pip` instead of `pip3` depending on your Linux distribution):
```
pip3 install notebook
```

### pystan 

To install `pystan`, run the follwing shell command:
```
pip3 install pystan
```


## Installation instructions for Windows

We recommend to install Anaconda via the instructions [here](https://docs.anaconda.com/anaconda/install/). The Anaconda installation will provide Jupyter notebooks along with Python, including the `numpy` and `scipy` packages.

### pystan 

Once Anaconda is installed, open the Anaconda Prompt and type (see [here](https://docs.anaconda.com/anaconda/user-guide/tasks/install-packages/) for more information):
```
conda install pystan 
```


