## 0.10.18 (May 07, 2024)
  - [FIX]: only attempt to create the marginal photometry plot if photometry has been printed to the output file

## 0.10.17 (March 26, 2024)
  - [ADD] Added option to print LaTeX table in column-wise order (one parameter per row)
  - [MOD] Changing the default credible interval to be printed to the LaTeX table to 68%
  - [ADD] Added command line option to average or not the positive/negative errors printed in the output LaTeX table
  - [FIX]Fixing computationof LaTeX table when the summary catalogue has been computed flattening the columns
  - [FIX] Fixing command line options dependencies: some options require the  option

## 0.10.16 (March 14, 2024)
  - [FIX] Fixed y-axis label (units) when plotting EWs vs line fluxes
  - [FIX] Spectral indices labels were not correctly printed when reading them from a JSON file

## 0.10.15 (February 27, 2024)
  - [FIX] Fixing error when line labels in the spectral indices configuration file do not match the keys in the JSON file

## 0.10.14 (February 27, 2024)
  - [FEAT] Adding feature to correctly plot upper/lower limits as uppper/lower facing triangles
  - [ADD] Adding option to specify the definition of the upper/lower limits in units of sigma

## 0.10.13 (February 13, 2024)
  - [FIX] Missing default (False) flag for --compute-summary option

## 0.10.12 (January 26, 2024)
  - [FIX] Fixing errors when no summary_config.json file is present

## 0.10.11 (January 24, 2024)
  - [FIX] Fixing passing a summary_config.json file to configure the summary catalogue computation

## 0.10.10 (January 24, 2024)
  - [FIX] Fixing case column names matching

## 0.10.9 (December 06, 2023)
  - [FIX] Fixing errors of filters with names longer than 20 characters

## 0.10.8 (November 22, 2023)
  - [ADD] Reordering the hdu extensions so that the POSTERIOR PDF extension is always checked first
  - [FIX] Other fix for upper-lower cases mismatches. [ADD] Also added a sanity check that all parameters for which we require a LaTeX table have been printed
  - [FIX] Fixing column names mismatch caused by upper/lower cases
  - [ADD] Added horizontal line to show the level of 0 flux
  - [ADD] Added two command line options to set the rotation of the x-axis labels and to set the y-range of the residuals

## 0.10.7 (November 10, 2023)
  - [FIX] Fixing error `TypeError: Cannot compare structured or void to non-void arrays.` and removing parameters not present in the result file being analyzed

## 0.10.6 (November 10, 2023)
  - [ADD] Allowing the params_names.json file to contain entries not present in the output file, for instance adjustable parameters not fitted for. This enables users to use a common params_names.json file for different Beagle parameter files

## 0.10.5 (November 10, 2023)
  - Merge branch 'feature/log_level' into develop
  - [ADD] Adding possibility to set logging level from command line

## 0.10.4 (October 13, 2023)
  - Merge branch 'fix/crontab' into develop
  - [ADD] Added check so that when --show-calibration-option is set the presence of a  calibration polynomial is checked, to avoid encountering errors (this enable using the same pyp_beagle command for runs w and w/o calibration polynomial)
  - [FIX] numpyp.int --> int

## 0.10.3 (August 08, 2023)
  - Skip objects with all fluxes set to 0
  - Avoid using NaN or Inf as axes limits, fixes #54
  - Correctly dealing with observed spectra with negative errors, i.e. we now ignore the wl bins corresponding to negative errors

## 0.10.2 (April 28, 2023)
  - Removing undefined float64 data type

## 0.10.1 (April 25, 2023)
  - Version bump 0.10.0
  - Removed spurious print statements
  - Python 3.7 --> 3.8
  - Refactoring to the latest versions of numpy and matplotlib
  - Removing mandatory presence of json file to configure emission line labels
  - Enable units definition in spectral indices catalogue
  - Version bump 0.10.0
  - Python 3.7 --> 3.8
  - Refactoring to the latest versions of numpy and matplotlib
  - Removing mandatory presence of json file to configure emission line labels
  - Enable units definition in spectral indices catalogue
  - Version bump 0.9.13
  - Several cosmetic improvements to the way residuals are shown, in particular by introducing the possibility of showing residuals in absolute units, relative units, or in units of sigma
  - Show tick marks at integer locations
  - Added option to show different type of residuals: absolute, relative, or in units of sigma (default option)
  - Version bump 0.9.12
  - Bug fix: the '--show-line-labels' option was still creating issues for labels at the edges of the spectra
  - Version bump 0.9.11
  - Bug fix: emission lines labelling in marginal plots now works again
  - Updated requirements
  - Version bump 0.9.10
  - The options '--show-full-SED' and '--show-MAP-SED' now work well for the spectra marginal plots. A new '--n-full-SED' option has been added to define the number of full SED to overplot on the marginal plots
  - Version bump 0.9.9
  - Bug fix: residuals can be correctly computed even when the marginal wl array and the data wl array do not have the same size, fixes #52
  - Version bump 0.9.8
  - Removing unsupported 'papertype' option from matplotlib fig.savefig method
  - Update README.md
  - Version bump 0.9.7
  - Bug fix: an issue was introduced when adding support for the MAP solution in the summary catalogue, this fix closes #51
  - Bug fix: correctly dealing with log flux scale in the presence of negative or zero fluxes

## 0.10.0 (April 25, 2023)
  - Removed spurious print statements
  - Python 3.7 --> 3.8
  - Refactoring to the latest versions of numpy and matplotlib
  - Removing mandatory presence of json file to configure emission line labels
  - Enable units definition in spectral indices catalogue

## 0.9.13 (November 16, 2022)
  - Several cosmetic improvements to the way residuals are shown, in particular by introducing the possibility of showing residuals in absolute units, relative units, or in units of sigma
  - Show tick marks at integer locations
  - Added option to show different type of residuals: absolute, relative, or in units of sigma (default option)

## 0.9.12 (November 15, 2022)
  - Bug fix: the '--show-line-labels' option was still creating issues for labels at the edges of the spectra

## 0.9.11 (November 15, 2022)
  - Bug fix: emission lines labelling in marginal plots now works again
  - Updated requirements

## 0.9.10 (November 13, 2022)
  - The options '--show-full-SED' and '--show-MAP-SED' now work well for the spectra marginal plots. A new '--n-full-SED' option has been added to define the number of full SED to overplot on the marginal plots

## 0.9.9 (November 11, 2022)
  - Bug fix: residuals can be correctly computed even when the marginal wl array and the data wl array do not have the same size, fixes #52

## 0.9.8 (November 11, 2022)
  - Removing unsupported 'papertype' option from matplotlib fig.savefig method
  - Update README.md

## 0.9.7 (November 11, 2022)
  - Bug fix: an issue was introduced when adding support for the MAP solution in the summary catalogue, this fix closes #51
  - Bug fix: correctly dealing with log flux scale in the presence of negative or zero fluxes

## 0.9.6 (July 22, 2022)
  - Raise error whne unable to match object IDs and input spectra
  - Allow post-processing of spectra with redshift as an adjustable parametter, i.e. the redshift column can contain different values
  - Defining backend when on Mac OS

## 0.9.5 (July 22, 2022)
  - Adding some logging messages to ease debugging

## 0.9.4 (July 21, 2022)
  - Added new fields to the summary catalogue: the MAP value for each parameter, and the probability, likelihood and chi_square corresponding to the MAP solution. Also added the number of data points used in the fitting to ease the calculation of a pseudo reduceed chi_square. Closes #49
  - Added command-line options to overwrite existing files and to flatten the columns in the output FITS files
  - Replacing deprecated Astropy FITS clobber keyword, fixes #50

## 0.9.3 (May 03, 2022)
  - Fixing plotting of full SED (spectrum) on marginal plots for photometry. Units conversion (F_lambda to F_nu) was wrong (a regression that appeared at some point). Fixes #47
  - Remove Matplotlib deprecated option (see https://matplotlib.org/3.3.1/gallery/text_labels_and_annotations/usetex_baseline_test.html), fixes #48

## 0.9.2 (October 29, 2021)
  - Enabling setting wl units for photometry plots as well (via the --wl-units option)

## 0.9.1 (June 28, 2021)
  - Bug fix: now supporting negative fluxes (e.g. absorption lines)
  - Removed dependency on WeightedKDE since scipy.stats now supports weighted gaussian_kde

## 0.9.0 (June 14, 2021)
  - Merge pull request #46 from eclake/calibration_branch_new
  - removing print statements and fixing bugs
  - Update beagle_spectra.py
  - Update beagle_spectra.py
  - Update beagle_calibration_correction.py
  - adding calibration_correction plot

## 0.8.3 (May 04, 2021)
  - Removed dependencies which are already included in the Python Standard Library

## 0.8.2 (April 26, 2021)
  - Fixing wrong module imports, fixes #44

## 0.8.1 (April 26, 2021)
  - Modifying Github action badge
  - Editing Github workflow
  - Modifying the workflow to trigger the action at releases
  - Adding workflow badge
  - Changed Github workflow path and file name

## 0.8.0 (April 26, 2021)
  - Merge branch 'feature/github_actions' into develop
  - Removing Shippable file
  - Modifying Github actions workflow
  - Editing Github workflow/actions
  - Updated requirements
  - Adding directory for Github actions

## 0.7.4 (March 23, 2021)
  - Fixing some Python 3 conversion issues
  - Corrected a wrong import
  - Enabling reading the ID column name of a photometric catalogue from the filters configuration file (it can still be overridden through the command line argument)

## 0.7.3 (March 22, 2021)
  - Trying to solve CI error: can't find Rust compiler

## 0.7.2 (March 22, 2021)
  - RuntimeError: latex was not able to process the following string:b'lp' solved following https://stackoverflow.com/a/32020370, i.e. redefining the minus sign in the preamble 
  - Recognizing use of env variables in the filters definition; skipping the header line in the ASCII filter transmission files
  - Trying to solve a LaTeX error RuntimeError: latex was not able to process the following string: b'lp' using hints from here https://github.com/garrettj403/SciencePlots/issues/2 and here https://stackoverflow.com/a/53080504
  - Adding LaTex use through our matplotlib configuration
  - Updated env variables
  - Updating README
  - Adding missing import; adding configparser option since the default option has changed

## 0.7.1 (January 21, 2021)
  - Python 3.8 not yet available on Shippable, using 3.7 instead
  - Updating README to reflect move to Python 3

## 0.7.0 (January 21, 2021)
  - Merge branch 'python_3' into develop
  - Moving to support of Python 3 (3.8) exclusively
  - Updated setup for Python 3
  - Updated modules for Python 3

## 0.6.6 (July 16, 2019)
  - Fix pypi failure: specifying that the README is in Markdown format, otherwise pushing to pypi fails, see https://pypi.org/help/#description-content-type and https://github.com/di/markdown-description-example
  - Bug fix: the marginal photometry plot was failing when the filter configuration file contained bands not present in the observations, but only computed for the models
  - Added option to specify flux units for photometry
  - Update README.md

## 0.6.5 (April 05, 2019)
  - Read object redshift from the Beagle output file, so that pyp_beagle correctly deals with redshifts when provided in the prior catalogue; correctly plots masked regions read from the 'MARGINAL SED MASK' extension
  - Correctly identify object ID from Beagle output file name
  - Cosmetics improvements to the marginal photometry plot

## 0.6.4 (March 27, 2019)
  - Introduced option to overplot on a photometric SED the maximum-a-posteriori (MAP) spectrum
  - Exploiting the *args and **kwargs appraoches to simplify the passing of command-line arguments to the different classes

## 0.6.3 (March 08, 2019)
  - Bug fix: file names of spectra and object IDs were not correctly matched

## 0.6.2 (February 28, 2019)
  - Correctly dealing with lines not observed (i.e. negative fluxes and errors in the observed catalogue)

## 0.6.1 (February 28, 2019)
  - Implemented option to print the values of observed and model fluxes on the plot
  - Removed unused function
  - Changed parser variable name
  - Cosmetics: added dotted lines linking symbols in the plot to the line labels on the x-axis

## 0.6.0 (February 27, 2019)
  - Merge branch 'feature/spectral_indices' into develop
  - Added SpectralIndices module
  - Finished integration of module to plot spectral indices (integrated line fluxes and EWs)
  - Added option to print values of line fluxes in the output image
  - Started integrating the new 'SpectralIndices' class in the main module
  - Started implementing the new module to post-process fitting to spectral indices / integrated fluxes / EWs
  - Using the new class 'ObservedCatalogue'
  - Added new module containing a generic 'observed catalogue'
  - Changed mail address
  - Removed support for python 2.6

## 0.5.20 (October 02, 2018)
  - The font size in the triangle plots can now be set using the '--fontsize' command-line option; for aesthetical reasons, reduced the number of tick marks and tick labels in the triangle plots
  - In the absence of a JSON file, the summary catalogue is now computed by default for the 'POSTERIOR PDF' extension

## 0.5.19 (September 05, 2018)
  - Removing LaTeX package aas_macros.sty, deals with #34

## 0.5.18 (September 05, 2018)
  - Included matplotlib configuration directly inside the package, to avoid the need of having to customize the matplotlibrc configuration file

## 0.5.17 (August 28, 2018)
  - Specifying astropy version, otherwise the installation will fail, closes #33

## 0.5.16 (August 02, 2018)
  - Added option to show labels of photometric filters in the marginal plot

## 0.5.15 (August 01, 2018)
  - Bug fix: pyp_beagle was crashing when the filter throughputs curves were passed via the 'fileName:' token inside the filter configuration file, in the absence of the 'FILTERS THROUGHPUTS' FITS file

## 0.5.14 (July 30, 2018)
  - Updated dependencies
  - Allow spectra without a redshift header keyword to be correctly processed
  - Added version file

## 0.5.13 (July 06, 2018)
  - Introduced possibility to use multiprocessing to speed up the calculation of the summary catalogue; added natural sorting or rows in the summary catalogue
  - Added dependency on 'natsort' package for natural sorting
  - Renamed 'np' variable to 'n_proc' to avoid confusion with numpy ---> np

## 0.5.12 (June 04, 2018)
  - Enable use of IDs longer than 20 characters (100 characters is current limit, but it can be easily changed), deals with #32

## 0.5.11 (April 05, 2018)
  - Implemented option to plot log10(fluxes) for photometry

# Change Log

## [0.5.10](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.10) ("03/28/2018")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.9...0.5.10)

**Closed issues:**

- Meaning of keywords in the summary file [\#30](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/30)

## [0.5.9](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.9) ("12/13/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.8...0.5.9)

**Closed issues:**

- TypeError: 'AxesSubplot' object has no attribute '\_\_getitem\_\_' [\#28](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/28)

## [0.5.8](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.8) ("12/07/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.7...0.5.8)

## [0.5.7](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.7) ("11/29/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.6...0.5.7)

**Closed issues:**

- Error with 'gridspec\_kw' in version 0.5.6 [\#29](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/29)

## [0.5.6](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.6) ("09/06/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.5...0.5.6)

## [0.5.5](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.5) ("08/24/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.4...0.5.5)

**Closed issues:**

- Missing dependencies?  "No module named set\_shared\_labels" [\#27](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/27)

## [0.5.4](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.4) ("08/08/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.3...0.5.4)

## [0.5.3](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.3) ("08/08/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.2...0.5.3)

## [0.5.2](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.2) ("07/28/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.1...0.5.2)

## [0.5.1](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.1) ("07/27/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.5.0...0.5.1)

## [0.5.0](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.5.0) ("07/26/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.6...0.5.0)

**Closed issues:**

- Issues with --plot-triangle and --plot-marginal [\#25](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/25)

## [0.4.6](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.6) ("06/22/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.5...0.4.6)

## [0.4.5](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.5) ("06/22/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.4...0.4.5)

## [0.4.4](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.4) ("06/07/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.3...0.4.4)

## [0.4.3](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.3) ("05/16/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.2...0.4.3)

## [0.4.2](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.2) ("05/04/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.1...0.4.2)

## [0.4.1](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.1) ("05/02/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.4.0...0.4.1)

**Closed issues:**

- Error with --plot-marginal \(fit to multiple spectra\) [\#23](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/23)
- ImportError: No module named autoscale [\#22](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/22)
- Summary Catalog Error with "compute" [\#20](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/20)

## [0.4.0](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.4.0) ("03/30/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.3.1...0.4.0)

**Closed issues:**

- Error running ./postprocess\_beagle\_results.py ImportError: No module named WeightedKDE [\#21](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/21)
- BANGS\_summary\_catalogue masses [\#3](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/3)

## [0.3.1](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.3.1) ("03/29/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.3.0...0.3.1)

## [0.3.0](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.3.0) ("03/27/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.10...0.3.0)

## [0.2.10](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.10) ("03/08/2017")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.9...0.2.10)

**Closed issues:**

- numpy TypeError [\#19](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/19)

## [0.2.9](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.9) ("12/08/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.8...0.2.9)

## [0.2.8](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.8) ("11/01/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.7...0.2.8)

**Fixed bugs:**

- error creating triangle plot with single value parameter [\#17](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/17)

**Closed issues:**

- Error creating corner plots [\#16](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/16)
- new hotfix branch [\#15](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/15)
- running of the spectral example [\#11](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/11)

## [0.2.7](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.7) ("10/06/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.6...0.2.7)

**Closed issues:**

- Crash while making summary catalog [\#14](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/14)

## [0.2.6](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.6) ("10/05/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.5...0.2.6)

**Closed issues:**

- python packages for PyP-BEAGLE [\#7](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/7)

## [0.2.5](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.5) ("10/05/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.4...0.2.5)

**Closed issues:**

- matplotlibrc file [\#12](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/12)

## [0.2.4](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.4) ("10/05/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.3...0.2.4)

**Closed issues:**

- Crashin while making summary catalog [\#13](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/13)

## [0.2.3](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.3) ("09/01/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.2...0.2.3)

**Closed issues:**

- Error running the configuration summary [\#10](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/10)

## [0.2.2](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.2) ("08/20/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.1...0.2.2)

## [0.2.1](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.1) ("07/21/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.2.0...0.2.1)

## [0.2.0](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.2.0) ("07/21/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.7...0.2.0)

**Fixed bugs:**

- Calculations of medians and intervals in beagle\_summary\_catalogue [\#9](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/9)

**Closed issues:**

- latex - python incompatibility [\#8](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/8)

## [0.1.7](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.7) ("07/11/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.6...0.1.7)

## [0.1.6](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.6) ("05/13/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.5...0.1.6)

## [0.1.5](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.5) ("05/13/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.4...0.1.5)

**Closed issues:**

- Filter file parsing problem [\#5](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/5)
- Is there some way I can get an example for how to use PyP-BEAGLE? [\#4](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/4)

## [0.1.4](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.4) ("03/17/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.3...0.1.4)

## [0.1.3](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.3) ("03/16/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.2...0.1.3)

## [0.1.2](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.2) ("03/15/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.1...0.1.2)

## [0.1.1](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.1) ("03/11/2016")
[Full Changelog](https://github.com/jacopo-chevallard/PyP-BEAGLE/compare/0.1.0...0.1.1)

## [0.1.0](https://github.com/jacopo-chevallard/PyP-BEAGLE/tree/0.1.0) ("03/11/2016")
**Fixed bugs:**

- bangs\_photometry singular matrix [\#1](https://github.com/jacopo-chevallard/PyP-BEAGLE/issues/1)



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*
