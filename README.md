## Before using PyP-BEAGLE

* clone this repository!

* intall [Python Anaconda](https://www.continuum.io/downloads) (the Python 2.7 version)

* install the packages below, which are not included in the Anaconda Python distribution. You can install them with the command ``pip install <package name>``, where ``<package name>`` will be:

  * atpy

  * getdist

  * pathos

* copy the file ``script/matplotlibrc`` into your ``$HOME/.matplotlib`` folder (if the folder does not exist, create it). If you already have a customized  ``matplotlibrc`` file, then try to run the PyP-BEAGLE example files located in the ``tests`` folder, and use GNU ``diff`` command to find out which additions to your ``matplotlibrc`` are required to run PyP-BEAGLE.

* set the Beagle environment variables using the ``BEAGLE_env_variable.bash`` or ``BEAGLE_env_variable.csh`` located in the ``scripts`` folder. For this, you can simply add at the end of your ``.bashrc`` file the line
 ```bash
 source <path to the cloned PyP-BEAGLE repository>/scripts/BEAGLE_env_variable.bash
 ```

* set the ``PYP_BEAGLE`` environment variable to point to the location of the cloned ``PyP-BEAGLE`` repository

## Using PyP-BEAGLE

To start using PyP-BEAGLE, take a look at the example files in the tests folder. 

To run the ``test_photometry.py`` example on your output you will need to copy to the results folder the ``PyP-BEAGLE/files/params_names.json`` file, and edit it to reflect the ``fitted`` parameters set in your parameter file.

Running 
```shell
test_photometry.py --help
```
will then give you indications on how to proceed.
