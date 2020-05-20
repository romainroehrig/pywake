# pywake: A Python Interface to the wake parameterization

The present scripts implement an Python interface to the wake parameterization developped by Grandpeix and Lafore (2010), including more recent developments (e.g., prognostic wake density).

## Requirements

* gfortran 
* Python 2.7
* Python packages: netCDF4, numpy, matplotlib

## Installation and test

* `git clone https://github.com/romainroehrig/pywake`
* compile the wake library:
	* `cd wakelib`
	* `compile.sh` (you may need to modify this compile script)
* run the test: `python run_example_ARPEGE.py`
