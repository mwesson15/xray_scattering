# xray_scattering

# Introduction
This repository contains code to simulate X-Ray scattering off a sphere in ```sphere_3d.py``` and code to replicate the simulations of *Reconstructing detailed line profiles of lamellar gratings from GISAXS patterns with a Maxwell solver* (Soltwisch et. al. 2017) in ```soltwisch.py```.

These simulations are written in Python with [MEEP](https://meep.readthedocs.io/en/latest/), an open source finite-difference time-domain package for solving maxwell's equations.

# Installation
In order to run the code, first the environment needs to be set up.

We use [anaconda](https://docs.anaconda.com/anaconda/user-guide/faq/#anaconda-faq-35) with Python 3.6 as our base Python.
After installing anaconda, in order to create the environment run:
```
conda create -n meep -c conda-forge pymeep tqdm
```
This will create an environment named ```meep``` with the additional packages ```pymeep```(our MEEP installation) and ```tqdm```(a progress bar) installed. Prior to running the code activate the environment with:
```
source activate meep
```
Then you can run the code with:
```
python [filename]
```
When finished, you can deactivate the environment with:
```
source deactivate meep
```

In order to get the code, either clone this repo with:
```
git clone [repo-url]
```
Or copy and paste the code directly into your own file.

If you encounter any issues, the [MEEP installation instructions](https://meep.readthedocs.io/en/latest/) may be helpful.
