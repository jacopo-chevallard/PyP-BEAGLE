# Before using PyP-BEAGLE

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

# Using PyP-BEAGLE

## Spectroscopic data

The post-processing of Beagle results obtained by fitting spectroscopic data is performed by means of the script ``PyP-BEAGLE/postprocess_beagle_spectra.py``. You can visualize the possible command-line options of the script with the command
```csh
./postprocess_beagle_spectra.py --help
```

### Plotting the posterior probability distributions (aka "triangle plots")

#### Command

```csh
./postprocess_beagle_spectra.py -r <your Beagle results folder> \
--plot-triangle \
[-np <number of processors>] \
[--json-triangle <JSON triangle file>] \
[--mock-catalogue <input mock catalogue>] \
[--json-mock <JSON mock file>]
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<number of processors>`` is an integer indicating how many processors can be used for the parallel execution of the script. This is particularly important when producing plots for large (> 1000) samples, as the creation of each individual plot can take several tens of seconds.
* ``<JSON triangle file>`` is a JSON file used for the configuration of the triangle plot (which parameters should be plotted, log scale, plot limits, ...), an example can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/master/files/params_names.json);
* ``<input mock catalogue>`` indicates a Beagle FITS file containing the input (i.e. "true") physical parameters used to construct the noiseless SEDs which have then been fitted with Beagle (after the noise addition, which must be performed **outside** Beagle). Note that in this case, a ``<JSON mock file>`` must be passed, since we must instruct PyP-BEAGLE where (in which FITS extension and column) to find the "true" parameters. An example of the ``<JSON mock file>`` to be used in this case can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/master/files/params_names_mock.json).


#### Output

The successful execution of the script will create a set of ``*_triangle.pdf`` files (one per object) in the ``<your Beagle results folder>/pyp-beagle/plot`` folder.


### Plotting the comparison of data and model observables (aka "marginal plots")

#### Command

```csh
./postprocess_beagle_spectra.py -r <your Beagle results folder> \
--plot-marginal \
[-np <number of processors>] 
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<number of processors>`` is an integer indicating how many processors can be used for the parallel execution of the script. This is particularly important when producing plots for large (> 1000) samples, as the creation of each individual plot can take several tens of seconds.

#### Output

The successful execution of the script will create a set of ``*_marginal_SED_spec.pdf`` files (one per object) in the ``<your Beagle results folder>/pyp-beagle/plot`` folder.


### Computing a summary catalogue

#### Command

```csh
./postprocess_beagle_spectra.py -r <your Beagle results folder> 
--compute-summary
[--json-summary <JSON summary file>]
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<JSON summary file>`` is a JSON file used for the configuration of the summary catalogue, specifying for which parameters the summary statistics (posterior mean and median, 68 and 95 % credible regions) should be computed. An example can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/master/files/summary_config.json).

#### Output

The successful execution of the script will create the file ``<your Beagle results folder>/pyp-beagle/data/BEAGLE_summary_catalogue.fits``.

### Plotting the comparison of input and retrieved parameters when fitting mock observations

#### Command

```csh
./postprocess_beagle_spectra.py -r <your Beagle results folder> 
--mock-catalogue <input mock catalogue> \
--json-mock <JSON mock file>
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<input mock catalogue>`` indicates a Beagle FITS file containing the input (i.e. "true") physical parameters used to construct the noiseless SEDs which have then been fitted with Beagle (after the noise addition, which must be performed **outside** Beagle). Note that in this case, a ``<JSON mock file>`` must be passed, since we must instruct PyP-BEAGLE where (in which FITS extension and column) to find the "true" parameters. An example of the ``<JSON mock file>`` to be used in this case can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/master/files/params_names_mock.json).

#### Output

The successful execution of the script will create the files ``<your Beagle results folder>/pyp-beagle/plot/BEAGLE_mock_retrieved_params_hist.pdf`` and ``<your Beagle results folder>/pyp-beagle/plot/BEAGLE_mock_retrieved_params.pdf``.

