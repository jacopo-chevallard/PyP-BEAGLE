[![Upload Python Package](https://github.com/jacopo-chevallard/PyP-BEAGLE/actions/workflows/python-publish.yml/badge.svg)](https://github.com/jacopo-chevallard/PyP-BEAGLE/actions/workflows/python-publish.yml)

# PyP-BEAGLE

PyP-BEAGLE (**Py**thon **P**ostprocessing of **BEAGLE**) is a Python package to postprocess the analyses performed with the galaxy SED modelling tool [Beagle](http://www.jacopochevallard.org/beagle/) (**B**ay**E**sian **A**nalysis of **G**a**L**axy s**E**ds). PyP-BEAGLE allows one to create different types of publication-quality plots, LaTeX tables, as well as several higher level "summary" catalogues.  

# Installing PyP-BEAGLE

* Make sure that you have a (science-ready!) installation of Python 3.x (starting from PyP-BEAGLE version 0.7.0), for instance [Anaconda](https://www.continuum.io/downloads)

* To install PyP-BEAGLE simply run
  ```
  pip install pyp_beagle
  ```

# Known issues

* On a Mac OS, multiprocessing only works with the ``Agg`` backend. Make sure that your ``~/.matplotlib/matplotlibrc`` file contains the line
  ```
  backend      : Agg
  ```

* If you encounter errors related to LaTeX, or if the visual appearance of the plots is not satisfying, you can copy the matplotlib configuration file ``script/matplotlibrc`` into your ``$HOME/.matplotlib/`` folder (if the folder does not exist, create it). If you already have a customized  ``matplotlibrc`` file, then you can use the GNU ``diff`` command to update it.

* PyP-BEAGLE assumes that the Beagle environment variables are correctly set on your machine. Note that while these are the same environment variables used by Docker-Beagle (see [here](https://github.com/jacopo-chevallard/BEAGLE-general/wiki#running-beagle)), they have to point to the actual folders on your machine, not to the "virtual" folder that Docker-Beagle uses. To correctly set the environment variables, you can use the ``scripts/BEAGLE_env_variable.bash`` or ``scripts/BEAGLE_env_variable.csh`` files. In practice, after modifying the file to reflect your Beagle folder tree, you can simply add at the end of your ``.bashrc`` (or ``.tcshrc``, or equivalent) the line
 ```bash
 source <full path to the file>/BEAGLE_env_variable.bash
 ```

# Using PyP-BEAGLE

The post-processing of Beagle results is performed by means of the command ``pyp_beagle``. Since PyP-BEAGLE is often updated, you can visualize the (entire) possible options via the PyP-BEAGLE ``help``, with the command
```csh
pyp_beagle --help
```

Below we report a few of some common PyP-BEAGLE use cases and related commands.

* [triangle plots](#plotting-the-posterior-probability-distributions-aka-triangle-plots)
* [marginal plots](#plotting-the-comparison-of-data-and-model-observables-aka-marginal-plots)
* [summary catalogue](#computing-a-summary-catalogue)
* [true vs retrieved parameters](#plotting-the-comparison-of-input-and-retrieved-parameters-when-fitting-mock-observations)

### Plotting the posterior probability distributions (aka "triangle plots")

#### Command

```csh
pyp_beagle -r <your Beagle results folder> \
--plot-triangle \
[-np <number of processors>] \
[--json-triangle <JSON triangle file>] \
[--mock-catalogue <input mock catalogue>] \
[--json-mock <JSON mock file>]
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<number of processors>`` is an integer indicating how many processors can be used for the parallel execution of the script. This is particularly important when producing plots for large (> 1000) samples, as the creation of each individual plot can take several tens of seconds.
* ``<JSON triangle file>`` is a JSON file used for the configuration of the triangle plot (which parameters should be plotted, log scale, plot limits, ...), an example can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/0996fd3c6b271e15452b7edee6627bc7fbc68675/PyP-BEAGLE/files/params_names.json);
* ``<input mock catalogue>`` indicates a Beagle FITS file containing the input (i.e. "true") physical parameters used to construct the noiseless SEDs which have then been fitted with Beagle (after the noise addition, which must be performed **outside** Beagle). Note that in this case, a ``<JSON mock file>`` must be passed, since we must instruct PyP-BEAGLE where (in which FITS extension and column) to find the "true" parameters. An example of the ``<JSON mock file>`` to be used in this case can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/0996fd3c6b271e15452b7edee6627bc7fbc68675/PyP-BEAGLE/files/params_names_mock.json).


#### Output

The successful execution of the script will create a set of ``*_triangle.pdf`` files (one per object) in the ``<your Beagle results folder>/pyp-beagle/plot`` folder.


### Plotting the comparison of data and model observables (aka "marginal plots")

#### Command

```csh
pyp_beagle -r <your Beagle results folder> \
--plot-marginal \
[-np <number of processors>] \
[--log-wavelength] \
[--plot-line-labels] \
[--spectral-resolution <resolution>] \
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<number of processors>`` is an integer indicating how many processors can be used for the parallel execution of the script. This is particularly important when producing plots for large (> 1000) samples, as the creation of each individual plot can take several tens of seconds;
* ``<resolution>`` is a float indicating the resolution of the spectra, and it is used to determine which emission line labels are printed on the plot.

#### Output

The successful execution of the script will create a set of ``*_marginal_SED_spec.pdf`` files (one per object) in the ``<your Beagle results folder>/pyp-beagle/plot`` folder.


### Computing a summary catalogue

#### Command

```csh
pyp_beagle -r <your Beagle results folder> 
--compute-summary
[--json-summary <JSON summary file>]
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<JSON summary file>`` is a JSON file used for the configuration of the summary catalogue, specifying for which parameters the summary statistics (posterior mean and median, 68 and 95 % credible regions) should be computed. An example can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/0996fd3c6b271e15452b7edee6627bc7fbc68675/PyP-BEAGLE/files/summary_config.json).

#### Output

The successful execution of the script will create the file ``<your Beagle results folder>/pyp-beagle/data/BEAGLE_summary_catalogue.fits``.

#### Description

In the ``POSTERIOR PDF`` extension we have added some quantities related to the MAP = Maximum-a-Posteriori solution, namely the probability (``MAP_probability``), log-likelihood (``MAP_ln_likelihood``), chi-square (``MAP_chi_square``), and number of data points used in the fitting (``MAP_n_data``). These quantities enable a quick "frequentist-like" check of the goodness-of-the-fit of the MAP solution.

The physical parameters corresponding to the MAP solution are indicated as ``<parameter_name>_MAP`` (e.g. ``mass_MAP``).

### Plotting the comparison of input and retrieved parameters when fitting mock observations

#### Command

```csh
pyp_beagle -r <your Beagle results folder> 
--mock-catalogue <input mock catalogue> \
--json-mock <JSON mock file>
```

where
* ``<your Beagle results folder>`` must be replaced by the full path to the Beagle output directory;
* ``<input mock catalogue>`` indicates a Beagle FITS file containing the input (i.e. "true") physical parameters used to construct the noiseless SEDs which have then been fitted with Beagle (after the noise addition, which must be performed **outside** Beagle). Note that in this case, a ``<JSON mock file>`` must be passed, since we must instruct PyP-BEAGLE where (in which FITS extension and column) to find the "true" parameters. An example of the ``<JSON mock file>`` to be used in this case can be found [here](https://github.com/jacopo-chevallard/PyP-BEAGLE/blob/0996fd3c6b271e15452b7edee6627bc7fbc68675/PyP-BEAGLE/files/params_names_mock.json).

#### Output

The successful execution of the script will create the files ``<your Beagle results folder>/pyp-beagle/plot/BEAGLE_mock_retrieved_params_hist.pdf`` and ``<your Beagle results folder>/pyp-beagle/plot/BEAGLE_mock_retrieved_params.pdf``.

