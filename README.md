# atca_processing
General scripts that are useful for processing ATCA data and making lightcurves 

measureflux_casa.py:

should be run WITHIN casa! Not with casapy (or python) uses casa modules and for some reason only casa seems to handle the .cl files

processing.py:

run from terminal with casapy processing.py [options], note: you must use the casapy version or if you have the modular version of python install so it can call things like casatasks then ok, but otherwise, just use a casapy
