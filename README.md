#PyP-BEAGLE

**Important**: before running PyP-BEAGLE you must copy the file ``script/matplotlibrc`` into your ``$HOME/.matplotlib`` folder (if the folder does not exist, create it). If you already have a customized  ``matplotlibrc`` file, then try to run the PyP-BEAGLE example files located in the ``tests`` folder, and use GNU ``diff`` command to find out which additions to your ``matplotlibrc`` are required to run PyP-BEAGLE.


Required packages, not included in Anaconda Python, which can be installed through ``pip instal <package name>``:

* atpy

* getdist


## Using PyP-BEAGLE

To start using PyP-BEAGLE, take a look at the example files in the tests folder. 

To run the ``test_photometry.py`` example on your output you will to copy to the results folder the ``PyP-BEAGLE/files/params_names.json`` file, and edit it to reflect the ``fitted`` parameters set in your parameter file.

Running 
```shell
test_photometry.py --help
```
will then give you indications on how to proceed.
